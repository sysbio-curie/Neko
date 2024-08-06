import pandas as pd

from ... import data as _data


def huri(dataset: str = 'HI-union') -> pd.DataFrame:

    df = (
        _data.huri(dataset, translated = True).
        drop(['source', 'target'], axis = 1).
        rename(columns = lambda x: x.split('_')[0]).
        drop_duplicates()
    )

    return df

