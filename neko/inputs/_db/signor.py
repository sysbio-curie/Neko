from io import BytesIO
from typing import Mapping
from pathlib import Path
import logging
import os
import time
import uuid

import requests

import pandas as pd

from neko._cache import cache_dir

logger = logging.getLogger(__name__)

SIGNOR_URL = "https://signor.uniroma2.it/releases/getLatestRelease.php"
SIGNOR_ENTITIES_URL = "https://signor.uniroma2.it/download_complexes.php"
SIGNOR_CACHE_SUBDIR = 'signor'
SIGNOR_CACHE_FILENAME = 'SIGNOR_Human.tsv'
SIGNOR_DIRECT_VALUES = {
    'yes': True,
    'true': True,
    't': True,
    '1': True,
    'no': False,
    'false': False,
    'f': False,
    '0': False,
}
SIGNOR_REQUIRED_COLUMNS = {
    'IDA',
    'IDB',
    'DIRECT',
    'EFFECT',
    'ANNOTATOR',
    'PMID',
    'SIGNOR_ID',
}
SIGNOR_ENTITY_DOWNLOADS = {
    'complexes': {
        'filename': 'SIGNOR_complexes.csv',
        'submit': 'Download complex data',
        'columns': {'SIGNOR ID', 'COMPLEX NAME', 'LIST OF ENTITIES'},
    },
    'protein_families': {
        'filename': 'SIGNOR_protein_families.csv',
        'submit': 'Download protein family data',
        'columns': {'SIGNOR ID', 'PROT. FAMILY NAME', 'LIST OF ENTITIES'},
    },
    'phenotypes': {
        'filename': 'SIGNOR_phenotypes.csv',
        'submit': 'Download phenotype data',
        'columns': {
            'SIGNOR ID',
            'PHENOTYPE NAME',
            'PHENOTYPE DESCRIPTION',
        },
    },
    'stimuli': {
        'filename': 'SIGNOR_stimuli.csv',
        'submit': 'Download stimulus data',
        'columns': {
            'SIGNOR ID',
            'STIMULUS NAME',
            'STIMULUS DESCRIPTION',
        },
    },
}


def _cache_dir() -> Path:
    """Return the managed cache directory for SIGNOR resources."""

    return cache_dir(SIGNOR_CACHE_SUBDIR)


def _database_cache_path() -> Path:
    """Return the cached SIGNOR interaction-table path."""

    return _cache_dir() / SIGNOR_CACHE_FILENAME


def _write_dataframe_cache(
    df: pd.DataFrame,
    path: Path,
    separator: str,
) -> None:
    """Atomically write one validated SIGNOR DataFrame to the cache."""

    path.parent.mkdir(parents=True, exist_ok=True)
    token = f'{os.getpid()}.{uuid.uuid4().hex}'
    temporary_path = path.with_name(f'.{path.name}.{token}.tmp')

    try:
        df.to_csv(temporary_path, sep=separator, index=False)
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)


def _cache_dataframe(
    df: pd.DataFrame,
    path: Path,
    separator: str,
) -> None:
    """Best-effort cache write which never blocks use of downloaded data."""

    try:
        _write_dataframe_cache(df, path, separator)
    except OSError as error:
        logger.warning('Could not write SIGNOR cache %s: %s', path, error)


def _parse_signor_response(content: bytes) -> pd.DataFrame:
    """Parse and validate a SIGNOR TSV response."""

    try:
        df = pd.read_csv(BytesIO(content), sep='\t')
    except (pd.errors.ParserError, UnicodeDecodeError) as error:
        raise ValueError('SIGNOR returned an unreadable TSV response.') from error

    missing = SIGNOR_REQUIRED_COLUMNS - set(df.columns)

    if df.empty or missing:
        detail = (
            'the response was empty'
            if df.empty
            else f'missing columns: {", ".join(sorted(missing))}'
        )
        raise ValueError(
            f'SIGNOR returned an invalid TSV response ({detail}).',
        )

    return df


def _parse_signor_directness(direct: pd.Series) -> pd.Series:
    """Normalize every supported representation of SIGNOR directness."""

    direct_values = direct.astype('string').str.strip().str.casefold()
    is_direct = direct_values.map(SIGNOR_DIRECT_VALUES)
    invalid_direct = direct_values.notna() & is_direct.isna()

    if invalid_direct.any():
        invalid = sorted(set(direct_values.loc[invalid_direct].tolist()))
        raise ValueError(
            'Unrecognized values in SIGNOR DIRECT column: '
            f'{", ".join(invalid)}.',
        )

    return is_direct.fillna(False).astype(bool)


def _parse_signor_entity_response(
    content: bytes,
    required_columns: set[str],
) -> pd.DataFrame:
    """Parse and validate one of SIGNOR's semicolon-separated dictionaries."""

    try:
        df = pd.read_csv(BytesIO(content), sep=';')
    except (pd.errors.ParserError, UnicodeDecodeError) as error:
        raise ValueError(
            'SIGNOR returned an unreadable entity dictionary.',
        ) from error

    missing = required_columns - set(df.columns)

    if df.empty or missing:
        detail = (
            'the response was empty'
            if df.empty
            else f'missing columns: {", ".join(sorted(missing))}'
        )
        raise ValueError(
            f'SIGNOR returned an invalid entity dictionary ({detail}).',
        )

    return df.drop_duplicates().reset_index(drop=True)


def _read_cached_signor_database() -> pd.DataFrame | None:
    """Return a validated cached interaction table, if one is available."""

    path = _database_cache_path()

    if not path.exists():
        return None

    try:
        cached = _parse_signor_response(path.read_bytes())
        _parse_signor_directness(cached['DIRECT'])

        return cached
    except (OSError, ValueError) as error:
        logger.warning('Ignoring invalid SIGNOR cache %s: %s', path, error)

        return None


def _read_cached_signor_entity_dictionaries(
) -> dict[str, pd.DataFrame] | None:
    """Return all validated cached entity dictionaries, if complete."""

    cache = _cache_dir()
    dictionaries = {}

    for entity_type, config in SIGNOR_ENTITY_DOWNLOADS.items():
        path = cache / config['filename']

        if not path.exists():
            return None

        try:
            dictionaries[entity_type] = _parse_signor_entity_response(
                path.read_bytes(),
                config['columns'],
            )
        except (OSError, ValueError) as error:
            logger.warning(
                'Ignoring invalid SIGNOR entity cache %s: %s',
                path,
                error,
            )

            return None

    return dictionaries


def _load_signor_database() -> pd.DataFrame:
    """Load SIGNOR interactions from NeKo's cache or download and cache them."""

    cached = _read_cached_signor_database()

    if cached is not None:
        return cached

    downloaded = download_signor_database()
    _parse_signor_directness(downloaded['DIRECT'])
    _cache_dataframe(downloaded, _database_cache_path(), '\t')

    return downloaded


def _load_signor_entity_dictionaries() -> dict[str, pd.DataFrame]:
    """Load SIGNOR dictionaries from NeKo's cache or populate the cache."""

    cached = _read_cached_signor_entity_dictionaries()

    if cached is not None:
        return cached

    downloaded = download_signor_entity_dictionaries()
    cache = _cache_dir()

    for entity_type, config in SIGNOR_ENTITY_DOWNLOADS.items():
        _cache_dataframe(
            downloaded[entity_type],
            cache / config['filename'],
            ';',
        )

    return downloaded


def download_signor_entity_dictionaries(
    save_dir: str | Path | None = None,
    timeout: float = 120.0,
    attempts: int = 3,
    backoff: float = 0.5,
) -> dict[str, pd.DataFrame]:
    """
    Download SIGNOR dictionaries for non-protein interaction entities.

    The returned mapping contains ``complexes``, ``protein_families``,
    ``phenotypes``, and ``stimuli`` DataFrames. If ``save_dir`` is provided,
    the original semicolon-separated CSV responses are saved there as well.
    """

    if attempts < 1:
        raise ValueError('attempts must be at least 1.')

    destination = Path(save_dir) if save_dir is not None else None

    if destination is not None:
        destination.mkdir(parents=True, exist_ok=True)

    dictionaries = {}

    for entity_type, config in SIGNOR_ENTITY_DOWNLOADS.items():
        last_error = None

        for attempt in range(attempts):
            try:
                response = requests.post(
                    SIGNOR_ENTITIES_URL,
                    files={'submit': (None, config['submit'])},
                    timeout=timeout,
                )
                response.raise_for_status()
                entity_df = _parse_signor_entity_response(
                    response.content,
                    config['columns'],
                )
                dictionaries[entity_type] = entity_df

                if destination is not None:
                    (destination / config['filename']).write_bytes(
                        response.content,
                    )

                break

            except (requests.RequestException, ValueError) as error:
                last_error = error

                if attempt + 1 < attempts:
                    time.sleep(backoff * 2 ** attempt)
        else:
            raise RuntimeError(
                'Could not download a valid SIGNOR '
                f'{entity_type} dictionary after {attempts} attempts.',
            ) from last_error

    return dictionaries


def _dictionary_values(
    df: pd.DataFrame,
    value_column: str,
) -> dict[str, str]:
    """Build a clean SIGNOR-ID-to-value mapping from a dictionary table."""

    values = df[['SIGNOR ID', value_column]].dropna()

    return {
        str(identifier).strip(): str(value).strip()
        for identifier, value in values.itertuples(index=False, name=None)
    }


def _split_members(value: str) -> list[str]:
    """Split and clean a SIGNOR complex or family membership field."""

    return [member.strip() for member in str(value).split(',') if member.strip()]


def _signor_entity_mapping(
    dictionaries: Mapping[str, pd.DataFrame],
) -> dict[str, str]:
    """Create normalized identifiers for every SIGNOR-specific entity."""

    missing = set(SIGNOR_ENTITY_DOWNLOADS) - set(dictionaries)

    if missing:
        raise ValueError(
            'Missing SIGNOR entity dictionaries: '
            f'{", ".join(sorted(missing))}.',
        )

    complex_members = _dictionary_values(
        dictionaries['complexes'],
        'LIST OF ENTITIES',
    )
    family_members = _dictionary_values(
        dictionaries['protein_families'],
        'LIST OF ENTITIES',
    )

    def expand(identifier: str, ancestry: frozenset[str] = frozenset()):
        if identifier in ancestry:
            raise ValueError(
                f'Circular SIGNOR entity membership involving {identifier}.',
            )

        membership = complex_members.get(identifier)

        if membership is None:
            membership = family_members.get(identifier)

        if membership is None:
            return [identifier]

        expanded = []

        for member in _split_members(membership):
            expanded.extend(expand(member, ancestry | {identifier}))

        return expanded

    def unique_members(identifier: str) -> list[str]:
        return list(dict.fromkeys(expand(identifier)))

    mapping = {
        identifier: f'COMPLEX:{"_".join(unique_members(identifier))}'
        for identifier in complex_members
    }

    typed_names = (
        ('protein_families', 'PROT. FAMILY NAME', 'PROTEIN_FAMILY'),
        ('phenotypes', 'PHENOTYPE NAME', 'PHENOTYPE'),
        ('stimuli', 'STIMULUS NAME', 'STIMULUS'),
    )

    for dictionary_name, name_column, prefix in typed_names:
        names = _dictionary_values(
            dictionaries[dictionary_name],
            name_column,
        )
        mapping.update(
            {
                identifier: f'{prefix}:{name}'
                for identifier, name in names.items()
            },
        )

    return mapping


def normalize_signor_entities(
    df: pd.DataFrame,
    dictionaries: Mapping[str, pd.DataFrame],
) -> pd.DataFrame:
    """Replace proprietary SIGNOR endpoint IDs with normalized identifiers."""

    entity_mapping = _signor_entity_mapping(dictionaries)
    normalized = df.copy()

    for column in ('IDA', 'IDB'):
        normalized[column] = normalized[column].replace(entity_mapping)

    # The interaction export can briefly lead the downloadable dictionaries
    # when SIGNOR adds a new entity. In that case use the human-readable name
    # carried by the interaction itself while keeping a typed identifier. An
    # unresolved complex is deliberately called COMPLEX_NAME rather than
    # COMPLEX because its members are not available for OmniPath-style
    # expansion yet.
    fallback_prefixes = {
        'SIGNOR-C': 'COMPLEX_NAME',
        'SIGNOR-PF': 'PROTEIN_FAMILY',
        'SIGNOR-PH': 'PHENOTYPE',
        'SIGNOR-ST': 'STIMULUS',
    }

    for id_column, name_column in (
        ('IDA', 'ENTITYA'),
        ('IDB', 'ENTITYB'),
    ):
        if name_column not in normalized.columns:
            continue

        identifiers = normalized[id_column].astype('string')

        for signor_prefix, normalized_prefix in fallback_prefixes.items():
            mask = identifiers.str.startswith(signor_prefix, na=False)
            names = normalized.loc[mask, name_column].astype('string').str.strip()
            usable = names.notna() & names.ne('')
            indices = names.index[usable]
            normalized.loc[indices, id_column] = (
                normalized_prefix + ':' + names.loc[indices]
            )

    return normalized


def download_signor_database(
    save_path: str | None = None,
    timeout: float = 30.0,
    attempts: int = 3,
    backoff: float = 0.5,
) -> pd.DataFrame | None:
    """
    Download the latest stable SIGNOR human release.
    If save_path is provided, saves the file as TSV. Otherwise, returns a
    DataFrame.
    """
    if attempts < 1:
        raise ValueError('attempts must be at least 1.')

    last_error = None

    for attempt in range(attempts):

        try:
            response = requests.get(SIGNOR_URL, timeout=timeout)
            response.raise_for_status()
            df = _parse_signor_response(response.content)

            if save_path:
                with open(save_path, 'wb') as file:
                    file.write(response.content)

                return None

            return df

        except (requests.RequestException, ValueError) as error:
            last_error = error

            if attempt + 1 < attempts:
                time.sleep(backoff * 2 ** attempt)

    raise RuntimeError(
        f'Could not download a valid SIGNOR dataset after {attempts} attempts.',
    ) from last_error


def signor(
    path: str | None = None,
    entity_dictionaries: Mapping[str, pd.DataFrame] | None = None,
    normalize_entities: bool = True,
) -> pd.DataFrame:
    """
    SIGNOR database from a TSV, NeKo's cache, or a download.

    Processes a supplied local TSV. Without a path, uses NeKo's validated
    SIGNOR cache and downloads only resources missing from that cache.

    Parameters:
        path (str, optional):
            The path to the SIGNOR TSV. If None, uses the managed cache before
            downloading the database.
        entity_dictionaries (mapping, optional):
            Preloaded SIGNOR entity dictionaries. If omitted, the four
            dictionaries are loaded from the managed cache before downloading.
        normalize_entities (bool):
            Replace SIGNOR-specific complex, family, phenotype, and stimulus
            IDs with normalized, human-readable identifiers.

    Returns:
        pd.DataFrame: Processed SIGNOR interactions.
    """
    if path is None:
        df = _load_signor_database()
    else:
        df = pd.read_table(path)

    if normalize_entities:
        dictionaries = (
            entity_dictionaries
            if entity_dictionaries is not None
            else _load_signor_entity_dictionaries()
        )
        df = normalize_signor_entities(df, dictionaries)

    effects = df['EFFECT'].astype('string')
    unknown_effect = effects.str.strip().str.casefold().eq('unknown')
    df = df.loc[~unknown_effect.fillna(False)].copy()
    effects = df['EFFECT'].astype('string')

    is_direct = _parse_signor_directness(df['DIRECT'])

    # Transform the original dataframe into the desired format
    df = pd.DataFrame({
        'source': df['IDA'],
        'target': df['IDB'],
        # SIGNOR relationships are causal source-to-target relations. Its
        # DIRECT field describes whether the supporting evidence is direct,
        # rather than whether the graph edge has a direction.
        'is_directed': True,
        'is_direct': is_direct,
        'is_stimulation': effects.str.contains(
            'up-regulates', regex=False, na=False,
        ),
        'is_inhibition': effects.str.contains(
            'down-regulates', regex=False, na=False,
        ),
        'form_complex': effects.str.contains('complex', regex=False, na=False),
        'consensus_direction': True,
        'curation_effort': df['ANNOTATOR'],
        'references': df['PMID'],
        'sources': df['SIGNOR_ID'],
    })

    # Add the transformed DataFrame to the existing database
    df = _group_by_source_target(df)

    return df


def _group_by_source_target(df_ungrouped: pd.DataFrame) -> pd.DataFrame:

    if df_ungrouped.empty:
        result = df_ungrouped.copy()
        result['consensus_stimulation'] = pd.Series(dtype=bool)
        result['consensus_inhibition'] = pd.Series(dtype=bool)
        return result.reset_index(drop=True)

    group_columns = ['source', 'target']
    grouped = df_ungrouped.groupby(group_columns, sort=False, dropna=False)

    def join_strings(series):
        values = dict.fromkeys(
            str(value).strip()
            for value in series.dropna()
            if str(value).strip()
        )
        return '; '.join(values)

    votes = grouped[['is_stimulation', 'is_inhibition']].sum()
    result = grouped.aggregate({
        'is_directed': 'all',
        'is_direct': 'any',
        'is_stimulation': 'any',
        'is_inhibition': 'any',
        'form_complex': 'any',
        'consensus_direction': 'all',
        'curation_effort': join_strings,
        'references': join_strings,
        'sources': join_strings,
    }).reset_index()

    stimulation_votes = votes['is_stimulation'].to_numpy()
    inhibition_votes = votes['is_inhibition'].to_numpy()
    result['consensus_stimulation'] = stimulation_votes > inhibition_votes
    result['consensus_inhibition'] = inhibition_votes > stimulation_votes

    return result
