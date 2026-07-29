"""
Gene symbol <-> UniProt accession translation.

This module replaces the previous dependency on ``pypath.utils.mapping`` for
identifier translation. It combines two strategies:

1. A lazily-downloaded, disk-cached offline table covering all reviewed
   (Swiss-Prot) human UniProt entries, built from many small, bounded,
   retried requests to UniProt's paginated ``/uniprotkb/search`` endpoint.
   Once cached, lookups are local dict lookups with no network dependency.
2. A live fallback using the official job-based UniProt ID mapping REST API
   directly, for anything not found offline (non-human organisms, entries
   added after the cache was built, etc).

Neither strategy uses UniProt's ``/uniprotkb/stream`` endpoint, which streams
a single, long-lived, potentially huge response and has proven unreliable in
practice: a stalled connection there corrupts
the gzip stream being read incrementally, turning a transient network hiccup
into a confusing ``zlib`` error. The approaches used here issue many small,
independently-retried requests instead, so a single stalled request can be
retried in isolation without corrupting anything.

The cache is never populated at import time: it is only fetched the first
time a translation is actually requested.
"""

from __future__ import annotations

from pathlib import Path
import os
import re
import json
import time
import uuid
import hashlib
import logging

from requests.adapters import Retry, HTTPAdapter
import requests

from neko._cache import cache_dir

logger = logging.getLogger(__name__)

HUMAN_TAXON_ID = '9606'

_CACHE_SUBDIR = 'idmapping'
_TABLE_FILENAME = 'human_reviewed_idmapping.tsv'
_META_FILENAME = 'human_reviewed_idmapping.meta.json'
# Reviewed human proteome is ~20k entries; anything much lower indicates a
# truncated/incomplete download that should not be trusted.
_MIN_EXPECTED_ROWS = 10000

_SEARCH_URL = 'https://rest.uniprot.org/uniprotkb/search'
_RE_NEXT_LINK = re.compile(r'<([^>]+)>;\s*rel="next"')
_RE_UNIPROT_ACCESSION = re.compile(
    r'^[OPQ][0-9][A-Z0-9]{3}[0-9](-\d+)?$'
    r'|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(-\d+)?$'
)

_state = {
    'attempted': False,
    'symbol_to_uniprot': {},
    'uniprot_to_symbol': {},
    # In-memory memoization of live fallback lookups (including misses, as
    # `None`), so the same identifier is never re-queried over the network
    # more than once per process.
    'symbol_fallback_cache': {},
    'uniprot_fallback_cache': {},
}


def looks_like_uniprot_accession(value: str) -> bool:
    """Return True if `value` matches the UniProt accession pattern."""

    return bool(value) and bool(_RE_UNIPROT_ACCESSION.match(value))


def _cache_dir() -> Path:

    return cache_dir(_CACHE_SUBDIR)


def _table_path() -> Path:

    return _cache_dir() / _TABLE_FILENAME


def _meta_path() -> Path:

    return _cache_dir() / _META_FILENAME


def _new_session() -> requests.Session:

    session = requests.Session()
    retries = Retry(
        total = 5,
        backoff_factor = 0.5,
        status_forcelist = [429, 500, 502, 503, 504],
        allowed_methods = frozenset({'GET', 'POST'}),
    )
    session.mount('https://', HTTPAdapter(max_retries = retries))

    return session


def _get_next_link(headers) -> str | None:

    link = headers.get('Link')

    if link:
        match = _RE_NEXT_LINK.search(link)

        if match:
            return match.group(1)

    return None


def _fetch_reviewed_human_table(timeout: int = 30) -> list[tuple[str, str]]:
    """
    Download the accession / primary gene name table for all reviewed
    (Swiss-Prot) human UniProt entries, via small paginated requests.
    """

    session = _new_session()
    params = {
        'query': f'organism_id:{HUMAN_TAXON_ID} AND reviewed:true',
        'fields': 'accession,gene_primary',
        'format': 'tsv',
        'size': 500,
    }

    rows: list[tuple[str, str]] = []
    response = session.get(_SEARCH_URL, params = params, timeout = timeout)
    response.raise_for_status()

    def _parse(text: str) -> list[tuple[str, str]]:

        lines = text.splitlines()[1:]  # skip header
        parsed = []

        for line in lines:
            parts = line.split('\t')

            if len(parts) == 2 and parts[0]:
                parsed.append((parts[0], parts[1]))

        return parsed

    rows.extend(_parse(response.text))
    next_url = _get_next_link(response.headers)

    while next_url:
        response = session.get(next_url, timeout = timeout)
        response.raise_for_status()
        rows.extend(_parse(response.text))
        next_url = _get_next_link(response.headers)

    return rows


def _write_cache(rows: list[tuple[str, str]]) -> None:

    _validate_rows(rows)

    cache_dir = _cache_dir()
    cache_dir.mkdir(parents = True, exist_ok = True)

    table_path = _table_path()
    meta_path = _meta_path()
    token = f'{os.getpid()}.{uuid.uuid4().hex}'
    table_tmp_path = table_path.with_name(f'.{table_path.name}.{token}.tmp')
    meta_tmp_path = meta_path.with_name(f'.{meta_path.name}.{token}.tmp')
    table_text = 'Entry\tGene Names (primary)\n' + ''.join(
        f'{accession}\t{symbol}\n'
        for accession, symbol in rows
    )
    table_digest = hashlib.sha256(table_text.encode('utf-8')).hexdigest()

    try:
        with open(table_tmp_path, 'w', encoding = 'utf-8') as fh:
            fh.write(table_text)

        with open(meta_tmp_path, 'w', encoding = 'utf-8') as fh:
            json.dump(
                {
                    'rows': len(rows),
                    'fetched_at': time.time(),
                    'sha256': table_digest,
                },
                fh,
            )

        table_tmp_path.replace(table_path)
        meta_tmp_path.replace(meta_path)

    finally:
        table_tmp_path.unlink(missing_ok = True)
        meta_tmp_path.unlink(missing_ok = True)


def _validate_rows(rows: list[tuple[str, str]]) -> None:
    """Reject downloads too small to be a complete reviewed-human table."""

    if len(rows) < _MIN_EXPECTED_ROWS:
        raise ValueError(
            'Identifier mapping table is unexpectedly small '
            f'({len(rows)} rows; expected at least {_MIN_EXPECTED_ROWS}).',
        )


def _read_cache() -> list[tuple[str, str]] | None:
    """
    Read and validate the on-disk cache. Returns None if missing, unreadable,
    or failing basic sanity checks (so a truncated/corrupt cache is never
    silently trusted).
    """

    table_path = _table_path()
    meta_path = _meta_path()

    if not table_path.exists() or not meta_path.exists():
        return None

    try:
        with open(meta_path, encoding = 'utf-8') as fh:
            meta = json.load(fh)

        if not isinstance(meta, dict):
            raise ValueError('Identifier mapping metadata must be an object.')

        with open(table_path, encoding = 'utf-8') as fh:
            table_text = fh.read()

        lines = table_text.splitlines()

    except (OSError, ValueError) as e:
        logger.warning('Could not read identifier mapping cache: %s', e)

        return None

    rows = []

    for line in lines[1:]:
        parts = line.split('\t')

        if len(parts) == 2:
            rows.append((parts[0], parts[1]))

    digest = hashlib.sha256(table_text.encode('utf-8')).hexdigest()
    digest_mismatch = bool(meta.get('sha256')) and meta['sha256'] != digest

    if (
        len(rows) != meta.get('rows')
        or len(rows) < _MIN_EXPECTED_ROWS
        or digest_mismatch
    ):
        logger.warning(
            'Identifier mapping cache failed validation (%d rows); '
            'ignoring cached file.',
            len(rows),
        )

        return None

    return rows


def _build_dicts(
        rows: list[tuple[str, str]],
    ) -> tuple[dict[str, str], dict[str, str]]:

    symbol_to_uniprot: dict[str, str] = {}
    uniprot_to_symbol: dict[str, str] = {}

    for accession, symbol in rows:

        uniprot_to_symbol.setdefault(accession, symbol)

        if symbol:
            symbol_to_uniprot.setdefault(symbol, accession)

    return symbol_to_uniprot, uniprot_to_symbol


def _ensure_loaded(force_refresh: bool = False) -> None:

    if _state['attempted'] and not force_refresh:
        return

    _state['attempted'] = True

    rows = None if force_refresh else _read_cache()

    if rows is None:

        try:
            rows = _fetch_reviewed_human_table()
            _validate_rows(rows)
            _write_cache(rows)
        except Exception as e:
            logger.warning(
                'Could not download identifier mapping table (%s). '
                'Falling back to per-identifier live translation only.', e,
            )
            cached_rows = _read_cache() if force_refresh else None

            if cached_rows is not None:
                rows = cached_rows
            elif force_refresh and (
                _state['symbol_to_uniprot'] or _state['uniprot_to_symbol']
            ):
                # A failed explicit refresh must not discard a usable mapping
                # that was already loaded in this process.
                return
            else:
                rows = []

    symbol_to_uniprot, uniprot_to_symbol = _build_dicts(rows)
    _state['symbol_to_uniprot'] = symbol_to_uniprot
    _state['uniprot_to_symbol'] = uniprot_to_symbol


def refresh_cache(force: bool = True) -> None:
    """Force a fresh download of the offline identifier mapping table."""

    _ensure_loaded(force_refresh = force)


def _fallback_translate(
        ids: set[str],
        source: str,
        dest: str,
        taxon_id: str | None,
        poll_interval: float = 1.0,
        timeout: float = 20.0,
        request_timeout: float = 10.0,
    ) -> dict[str, str]:
    """
    Live single-shot translation using UniProt's job-based ID mapping REST
    API (submit / poll status / fetch results), for identifiers not found
    offline.

    This intentionally talks to the REST API directly with `requests` (and
    an explicit timeout on every single request) rather than going through
    `unipressed`'s convenience wrappers, whose internal `requests.get()`
    calls do not set a timeout at all: on this project's network, a single
    stalled connection there would block forever instead of failing fast
    and letting the caller move on.
    """

    session = _new_session()

    try:
        response = session.post(
            'https://rest.uniprot.org/idmapping/run',
            data = {
                'ids': ','.join(ids),
                'from': source,
                'to': dest,
                **({'taxId': taxon_id} if taxon_id else {}),
            },
            timeout = request_timeout,
        )
        response.raise_for_status()
        job_id = response.json()['jobId']

        deadline = time.monotonic() + timeout

        while time.monotonic() < deadline:
            status_response = session.get(
                f'https://rest.uniprot.org/idmapping/status/{job_id}',
                timeout = request_timeout,
            )
            status_response.raise_for_status()
            status_data = status_response.json()

            if 'results' in status_data or 'failedIds' in status_data:
                break

            job_status = status_data.get('jobStatus')

            if job_status == 'ERROR':
                return {}

            time.sleep(poll_interval)
        else:
            return {}

        results_response = session.get(
            f'https://rest.uniprot.org/idmapping/results/{job_id}',
            timeout = request_timeout,
        )
        results_response.raise_for_status()
        results = results_response.json().get('results', [])

        return {item['from']: item['to'] for item in results}

    except Exception as e:
        logger.warning('Live identifier translation failed: %s', e)

        return {}


def to_uniprot(label: str, organism: str = HUMAN_TAXON_ID) -> str | None:
    """
    Translate a gene symbol to its primary UniProt accession.

    If `label` is not recognized as a gene symbol but already looks like a
    valid UniProt accession, it is returned unchanged. Returns None if no
    translation could be found.
    """

    if not label:
        return None

    _ensure_loaded()

    accession = _state['symbol_to_uniprot'].get(label)

    if accession:
        return accession

    # Already looks like a UniProt accession (not a gene symbol): no need to
    # ask the live API to translate it as if it were one, just echo it back.
    if looks_like_uniprot_accession(label):
        return label

    if label in _state['symbol_fallback_cache']:
        return _state['symbol_fallback_cache'][label]

    dest = 'UniProtKB-Swiss-Prot' if organism == HUMAN_TAXON_ID else 'UniProtKB'
    fallback = _fallback_translate({label}, 'Gene_Name', dest, organism)
    result = fallback.get(label)
    _state['symbol_fallback_cache'][label] = result

    return result


def to_genesymbol(uniprot_id: str) -> str | None:
    """
    Translate a UniProt accession to its primary gene symbol.

    If `uniprot_id` is not recognized as an accession but does not look like
    one either (e.g. it is already a gene symbol), it is returned unchanged.
    Returns None if no translation could be found.
    """

    if not uniprot_id:
        return None

    _ensure_loaded()

    symbol = _state['uniprot_to_symbol'].get(uniprot_id)

    if symbol:
        return symbol

    # Does not look like a UniProt accession (presumably already a gene
    # symbol): no need to ask the live API to translate it as if it were
    # one, just echo it back.
    if not looks_like_uniprot_accession(uniprot_id):
        return uniprot_id

    if uniprot_id in _state['uniprot_fallback_cache']:
        return _state['uniprot_fallback_cache'][uniprot_id]

    fallback = _fallback_translate(
        {uniprot_id}, 'UniProtKB_AC-ID', 'Gene_Name', None,
    )
    result = fallback.get(uniprot_id)
    _state['uniprot_fallback_cache'][uniprot_id] = result

    return result
