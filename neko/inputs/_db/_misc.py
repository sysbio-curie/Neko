from typing import TYPE_CHECKING
import pandas as pd


def bool_col(df: pd.DataFrame, col: str, mappings: dict = None) -> pd.DataFrame:

    mappings = mappings or {
        "True": True,
        "False": False,
        "true": True,
        "false": False,
        1: True,
        -1: False,
        0: False,
    }

    if col not in df.columns:
        df[col] = False

    # Avoid FutureWarning: Downcasting behavior in `replace` is deprecated
    df[col] = df[col].replace(mappings)
    df[col] = df[col].infer_objects(copy=False).astype(bool)

    return df


def split_effect(
        df: pd.DataFrame,
        col: str = 'effect',
        mappings: dict = None,
        inhibition: bool = False,
    ) -> pd.DataFrame:
    """
    Split effects into two column notation (is_stimulation and is_inhibition).

    Parameters:
        col:
            Name of the column with the effect information.
        mapping:
            Mapping for the values to boolean.
        inhibition:
            Negate the effect: do it where True means inhibition.

    Returns:
        The data frame with col removed and is_stimulation and is_inhibition
        columns added.
    """

    df = bool_col(df, col, mappings)

    if inhibition:
        df[col] = ~df[col]

    df['is_inhibition'] = ~df[col]
    df.rename({col: 'is_stimulation'}, inplace = True, axis = 1)

    return df


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
