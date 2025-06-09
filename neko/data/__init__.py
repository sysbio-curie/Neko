import urllib.parse
import requests
import pandas as pd


def phosphosite_kinase_substrate():
    url = 'https://www.phosphosite.org/downloads/Kinase_Substrate_Dataset.gz'
    try:
        df = pd.read_csv(url, compression='gzip', sep='\t', low_memory=False)
    except Exception:
        # fallback for direct download if pandas fails
        r = requests.get(url)
        r.raise_for_status()
        from io import BytesIO
        df = pd.read_csv(BytesIO(r.content), compression='gzip', sep='\t', low_memory=False)
    return df


def phosphosite_regulatory_sites():
    url = 'https://www.phosphosite.org/downloads/Regulatory_sites.gz'
    try:
        df = pd.read_csv(url, compression='gzip', sep='\t', low_memory=False)
    except Exception:
        r = requests.get(url)
        r.raise_for_status()
        from io import BytesIO
        df = pd.read_csv(BytesIO(r.content), compression='gzip', sep='\t', low_memory=False)
    return df


def huri(dataset: str = 'HI-union', translated=True) -> pd.DataFrame:
    url = (
        'https://github.com/sysbio-curie/Medulloblastoma_project/'
        'raw/main/Huri_analysis/data/%s%s.csv'
    )
    url = url % (dataset, '_translated' if translated else '')
    df = pd.read_csv(url)
    return df
