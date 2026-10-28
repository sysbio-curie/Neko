from typing import TYPE_CHECKING

import numpy as np
import pandas as pd


_TRUE_STRINGS = {'true', 't', 'yes', 'y', '1', '1.0'}
_FALSE_STRINGS = {'false', 'f', 'no', 'n', '0', '0.0', '-1', '-1.0'}

# Effect label -> (is_stimulation, is_inhibition, form_complex)
_EFFECT_LABELS = {
    'stimulation': (True, False, False),
    'stimulate': (True, False, False),
    'stimulates': (True, False, False),
    'activation': (True, False, False),
    'activate': (True, False, False),
    'activates': (True, False, False),
    'positive': (True, False, False),
    '+': (True, False, False),
    'inhibition': (False, True, False),
    'inhibit': (False, True, False),
    'inhibits': (False, True, False),
    'repression': (False, True, False),
    'repress': (False, True, False),
    'represses': (False, True, False),
    'negative': (False, True, False),
    '-': (False, True, False),
    'bimodal': (True, True, False),
    'both': (True, True, False),
    'undefined': (False, False, False),
    'unknown': (False, False, False),
    'form complex': (False, False, True),
    'complex': (False, False, True),
}


def _is_missing(value) -> bool:

    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip()
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _to_bool(value, mappings: dict | None = None):
    """Boolean for a flag value, or ``None`` if the value is not recognised."""

    if mappings:
        try:
            if value in mappings:
                return bool(mappings[value])
        except TypeError:
            pass
    if _is_missing(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return value > 0
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in _TRUE_STRINGS:
            return True
        if normalized in _FALSE_STRINGS:
            return False
    return None


def _raise_unrecognised(col: str, values: list, accepted: str) -> None:

    shown = ', '.join(repr(v) for v in dict.fromkeys(values))
    raise ValueError(
        f"Column '{col}' contains unrecognised value(s): {shown}. "
        f"Accepted values: {accepted}."
    )


def bool_col(df: pd.DataFrame, col: str, mappings: dict = None) -> pd.DataFrame:
    """
    Convert a column to booleans, adding it as all False if it is missing.

    Accepts booleans, numbers (positive means True), and the strings
    true/false, t/f, yes/no, y/n, 1/0/-1 in any case. Missing values are False.

    Raises:
        ValueError: If the column holds any other value, for example a
            label such as ``"inhibition"`` that would otherwise be read
            as True.
    """

    if col not in df.columns:
        df[col] = False
        return df

    converted = [_to_bool(value, mappings) for value in df[col]]
    unrecognised = [v for v, c in zip(df[col], converted) if c is None]
    if unrecognised:
        _raise_unrecognised(
            col, unrecognised, 'True/False, yes/no, t/f, y/n or 1/0/-1',
        )
    df[col] = pd.Series(converted, index=df.index, dtype=bool)

    return df


def _effect_flags(value, inhibition: bool):

    if _is_missing(value):
        return _EFFECT_LABELS['undefined']
    if isinstance(value, str):
        label = value.strip().lower().replace('_', ' ')
        if label in _EFFECT_LABELS:
            return _EFFECT_LABELS[label]
    stimulation = _to_bool(value)
    if stimulation is None:
        return None
    if inhibition:
        stimulation = not stimulation
    return (stimulation, not stimulation, False)


def split_effect(
        df: pd.DataFrame,
        col: str = 'effect',
        mappings: dict = None,
        inhibition: bool = False,
    ) -> pd.DataFrame:
    """
    Split effects into two column notation (is_stimulation and is_inhibition).

    Recognised labels (case-insensitive): stimulation/activation/+,
    inhibition/repression/-, bimodal/both, undefined/unknown and
    form complex. Boolean-like values (True/False, 1/-1, yes/no) mean
    stimulation when true and inhibition when false. Missing values are
    undefined.

    Parameters:
        col:
            Name of the column with the effect information.
        mapping:
            Mapping for the values to boolean, applied before the labels.
        inhibition:
            Negate boolean-like values: do it where True means inhibition.

    Returns:
        The data frame with col removed and is_stimulation and is_inhibition
        columns added. Complex-formation labels also set ``form_complex``.

    Raises:
        ValueError: If any effect value is not recognised.
    """

    values = df[col]
    if mappings:
        values = values.map(lambda v: mappings.get(v, v) if _hashable(v) else v)
    flags = [_effect_flags(value, inhibition) for value in values]
    unrecognised = [v for v, f in zip(df[col], flags) if f is None]
    if unrecognised:
        _raise_unrecognised(
            col,
            unrecognised,
            'stimulation/activation/+, inhibition/repression/-, '
            'bimodal, undefined, form complex, or True/False and 1/-1',
        )

    df = df.drop(columns=col)
    df['is_stimulation'] = [f[0] for f in flags]
    df['is_inhibition'] = [f[1] for f in flags]
    complex_mask = pd.Series([f[2] for f in flags], index=df.index)
    if complex_mask.any():
        df = bool_col(df, 'form_complex')
        df['form_complex'] = df['form_complex'] | complex_mask

    return df


def _hashable(value) -> bool:

    try:
        hash(value)
    except TypeError:
        return False
    return True


def undirected_to_mutual(
        df: pd.DataFrame,
        source_col: str = 'source',
        target_col: str = 'target',
    ):

    cols = {source_col: target_col, target_col: source_col}

    df = (
        pd.concat([
            df,
            df.copy().rename(columns = cols),
        ]).
        drop_duplicates(
            subset = [source_col, target_col],
            ignore_index = True,
        )
    )

    return df
