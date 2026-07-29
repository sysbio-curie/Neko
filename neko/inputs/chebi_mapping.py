"""Lazy, cache-backed ChEBI accession-to-name lookup."""

from __future__ import annotations

from pathlib import Path
import logging
import os
import re
import uuid

from requests.adapters import HTTPAdapter, Retry
import requests

import pandas as pd

from neko._cache import cache_dir

logger = logging.getLogger(__name__)

CHEBI_COMPOUNDS_URL = (
    'https://ftp.ebi.ac.uk/pub/databases/chebi/'
    'flat_files/compounds.tsv.gz'
)

_CACHE_SUBDIR = 'chebi'
_COMPOUNDS_FILENAME = 'compounds.tsv.gz'
_MAPPING_FILENAME = 'chebi_names.tsv'
_MIN_EXPECTED_ROWS = 100_000
_CHEBI_PATTERN = re.compile(r'^CHEBI[:_](\d+)$', re.IGNORECASE)

_state = {
    'loaded': False,
    'id_to_name': {},
    'missing': set(),
}


def normalize_identifier(value: str | None) -> str | None:
    """Return a canonical ``CHEBI:<number>`` accession, when recognized."""

    if not isinstance(value, str):
        return None

    match = _CHEBI_PATTERN.fullmatch(value.strip())

    return f'CHEBI:{match.group(1)}' if match else None


def is_chebi_identifier(value: str | None) -> bool:
    """Return whether ``value`` is a ChEBI accession."""

    return normalize_identifier(value) is not None


def identifiers_in_frame(df: pd.DataFrame) -> set[str]:
    """Collect canonical ChEBI accessions from resource endpoints."""

    identifiers = set()

    for column in ('source', 'target'):
        if column not in df.columns:
            continue

        identifiers.update(
            normalized
            for value in df[column].dropna().unique()
            if (normalized := normalize_identifier(value)) is not None
        )

    return identifiers


def _cache_dir() -> Path:

    return cache_dir(_CACHE_SUBDIR)


def _compounds_path() -> Path:

    return _cache_dir() / _COMPOUNDS_FILENAME


def _mapping_path() -> Path:

    return _cache_dir() / _MAPPING_FILENAME


def _new_session() -> requests.Session:

    session = requests.Session()
    retries = Retry(
        total=3,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset({'GET'}),
    )
    session.mount('https://', HTTPAdapter(max_retries=retries))

    return session


def _read_mapping_cache() -> dict[str, str]:
    """Read the small derived mapping cache, ignoring invalid rows."""

    path = _mapping_path()

    if not path.exists():
        return {}

    try:
        df = pd.read_csv(path, sep='\t', dtype='string')
    except (OSError, pd.errors.ParserError, UnicodeDecodeError) as error:
        logger.warning('Could not read ChEBI name cache %s: %s', path, error)

        return {}

    required = ['chebi_accession', 'ascii_name']

    if set(required) - set(df.columns):
        logger.warning('Ignoring invalid ChEBI name cache %s.', path)

        return {}

    mapping = {}

    for identifier, name in df[required].itertuples(index=False, name=None):
        canonical = normalize_identifier(identifier)

        if canonical and pd.notna(name) and str(name).strip():
            mapping.setdefault(canonical, str(name).strip())

    return mapping


def _write_mapping_cache(mapping: dict[str, str]) -> None:
    """Atomically persist the derived ChEBI mapping."""

    path = _mapping_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    token = f'{os.getpid()}.{uuid.uuid4().hex}'
    temporary_path = path.with_name(f'.{path.name}.{token}.tmp')
    df = pd.DataFrame(
        sorted(mapping.items()),
        columns=['chebi_accession', 'ascii_name'],
    )

    try:
        df.to_csv(temporary_path, sep='\t', index=False)
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)


def _extract_names(path: Path, identifiers: set[str]) -> dict[str, str]:
    """Stream the compounds table and retain only requested accessions."""

    found = {}
    rows = 0
    wanted_columns = {'chebi_accession', 'ascii_name', 'name'}

    chunks = pd.read_csv(
        path,
        sep='\t',
        compression='gzip',
        usecols=lambda column: str(column).lower() in wanted_columns,
        chunksize=50_000,
        low_memory=False,
    )

    for chunk in chunks:
        rows += len(chunk)
        columns = {str(column).lower(): column for column in chunk.columns}
        identifier_column = columns.get('chebi_accession')
        ascii_column = columns.get('ascii_name')
        name_column = columns.get('name')

        if identifier_column is None or (
            ascii_column is None and name_column is None
        ):
            raise ValueError(
                'ChEBI compounds table is missing accession/name columns.',
            )

        matches = chunk[chunk[identifier_column].isin(identifiers)].copy()

        if matches.empty:
            continue

        labels = (
            matches[ascii_column].replace(r'^\s*$', pd.NA, regex=True)
            if ascii_column is not None
            else pd.Series(pd.NA, index=matches.index, dtype='string')
        )

        if name_column is not None:
            labels = labels.fillna(matches[name_column])

        for identifier, label in zip(matches[identifier_column], labels):
            canonical = normalize_identifier(identifier)

            if canonical and pd.notna(label) and str(label).strip():
                found.setdefault(canonical, str(label).strip())

    if rows < _MIN_EXPECTED_ROWS:
        raise ValueError(
            f'ChEBI compounds table is unexpectedly small ({rows} rows).',
        )

    return found


def _download_and_extract(
    identifiers: set[str],
    timeout: float = 120.0,
) -> dict[str, str]:
    """Download, validate, atomically cache, and query the compounds table."""

    path = _compounds_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    token = f'{os.getpid()}.{uuid.uuid4().hex}'
    temporary_path = path.with_name(f'.{path.name}.{token}.tmp')

    try:
        response = _new_session().get(CHEBI_COMPOUNDS_URL, timeout=timeout)
        response.raise_for_status()
        temporary_path.write_bytes(response.content)
        found = _extract_names(temporary_path, identifiers)
        temporary_path.replace(path)

        return found
    finally:
        temporary_path.unlink(missing_ok=True)


def _load_state() -> None:

    if not _state['loaded']:
        _state['id_to_name'] = _read_mapping_cache()
        _state['loaded'] = True


def _available_names(identifiers: set[str]) -> dict[str, str]:
    """Return the currently cached names for ``identifiers``."""

    return {
        identifier: _state['id_to_name'][identifier]
        for identifier in identifiers
        if identifier in _state['id_to_name']
    }


def ensure_names(identifiers: set[str]) -> dict[str, str]:
    """Ensure names for the requested ChEBI accessions are cached in memory."""

    canonical_ids = {
        canonical
        for value in identifiers
        if (canonical := normalize_identifier(value)) is not None
    }

    if not canonical_ids:
        return {}

    _load_state()
    missing = (
        canonical_ids
        - set(_state['id_to_name'])
        - set(_state['missing'])
    )

    if not missing:
        return _available_names(canonical_ids)

    path = _compounds_path()

    try:
        found = (
            _extract_names(path, missing)
            if path.exists()
            else _download_and_extract(missing)
        )
    except Exception as cached_error:
        if path.exists():
            logger.warning(
                'Ignoring invalid ChEBI compounds cache %s: %s',
                path,
                cached_error,
            )

            try:
                found = _download_and_extract(missing)
            except Exception as download_error:
                logger.warning('ChEBI name download failed: %s', download_error)
                _state['missing'].update(missing)

                return _available_names(canonical_ids)
        else:
            logger.warning('ChEBI name download failed: %s', cached_error)
            _state['missing'].update(missing)

            return _available_names(canonical_ids)

    _state['id_to_name'].update(found)
    _state['missing'].update(missing - set(found))

    if found:
        try:
            _write_mapping_cache(_state['id_to_name'])
        except OSError as error:
            logger.warning('Could not write ChEBI name cache: %s', error)

    return _available_names(canonical_ids)


def to_name(identifier: str) -> str | None:
    """Return the preferred ASCII name for a ChEBI accession, if available."""

    canonical = normalize_identifier(identifier)

    if canonical is None:
        return None

    ensure_names({canonical})

    return _state['id_to_name'].get(canonical)
