from typing import TYPE_CHECKING

if TYPE_CHECKING:

    import pandas as pd


def split_effect(
        df: pd.DataFrame,
        col: str,
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

    mappings = mapping or {
        "True": True,
        "False": False,
        "true": True,
        "false": False,
        1: True,
        -1: False,
        0: False,
    }
    df[col] = df[col].replace(mappings).astype(bool)

    if inhibition:
        df[col] = ~df[col]

    df['is_inhibition'] = ~df[col]
    df.rename({col: 'is_stimulation'}, inplace = True, axis = 1)

    return df
