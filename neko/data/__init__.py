import urllib.parse

import pandas as pd
from networkcommons.data.omics import _common as _ncdata


def _retrieve(path: str, ftype: str = 'tsv') -> pd.DataFrame:

    url = urllib.parse.urljoin(_ncdata._baseurl(), path)

    return _ncdata._open(url, ftype = ftype, df = True)


def phosphosite_kinase_substrate():

    return _retrieve('phosphosite/kinase-substrate.tsv')


def phosphosite_regulatory_sites():

    return _retrieve('phosphosite/regulatory-sites.tsv')


def huri(dataset: str = 'HI-union', translated = True) -> pd.DataFrame:

    url = (
        'https://github.com/sysbio-curie/Medulloblastoma_project/'
        'raw/main/Huri_analysis/data/%s%s.csv'
    )

    url = url % (dataset, '_translated' if translated else '')

    df = _ncdata._open(url, df = True)

    return df
