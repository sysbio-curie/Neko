import urllib.parse
import requests
import pandas as pd


def phosphosite_kinase_substrate():
    url = 'https://commons.omnipathdb.org/phosphosite/kinase-substrate.tsv'
    headers = {'User-Agent': 'Mozilla/5.0'}
    try:
        df = pd.read_csv(url, sep='\t', low_memory=False)
        return df
    except Exception:
        try:
            r = requests.get(url, headers=headers)
            r.raise_for_status()
            from io import StringIO
            df = pd.read_csv(StringIO(r.content.decode()), sep='\t', low_memory=False)
            return df
        except Exception as e:
            print("\n[NeKo] ERROR: Automatic download of PhosphoSitePlus Kinase_Substrate_Dataset failed.")
            print("Reason:", str(e))
            print("\nTo use this feature, please follow these steps:")
            print("1. Go to https://commons.omnipathdb.org/phosphosite/kinase-substrate.tsv\n2. Download 'kinase-substrate.tsv' manually.")
            print("3. Place the file in a directory of your choice (e.g., './neko/_data/').")
            print("4. When calling phosphosite(), pass the file path as the 'kinase_substrate' argument, e.g.:")
            print("   resources = phosphosite(kinase_substrate='./neko/_data/kinase-substrate.tsv')")
            raise


def phosphosite_regulatory_sites():
    url = 'https://commons.omnipathdb.org/phosphosite/regulatory-sites.tsv'
    headers = {'User-Agent': 'Mozilla/5.0'}
    try:
        df = pd.read_csv(url, sep='\t', low_memory=False)
        return df
    except Exception:
        try:
            r = requests.get(url, headers=headers)
            r.raise_for_status()
            from io import StringIO
            df = pd.read_csv(StringIO(r.content.decode()), sep='\t', low_memory=False)
            return df
        except Exception as e:
            print("\n[NeKo] ERROR: Automatic download of PhosphoSitePlus Regulatory_sites failed.")
            print("Reason:", str(e))
            print("\nTo use this feature, please follow these steps:")
            print("1. Go to https://commons.omnipathdb.org/phosphosite/regulatory-sites.tsv\n2. Download 'regulatory-sites.tsv' manually.")
            print("3. Place the file in a directory of your choice (e.g., './neko/_data/').")
            print("4. When calling phosphosite(), pass the file path as the 'regulatory_sites' argument, e.g.:")
            print("   resources = phosphosite(regulatory_sites='./neko/_data/regulatory-sites.tsv')")
            raise


def huri(dataset: str = 'HI-union', translated=True) -> pd.DataFrame:
    url = (
        'https://github.com/sysbio-curie/Medulloblastoma_project/'
        'raw/main/Huri_analysis/data/%s%s.csv'
    )
    url = url % (dataset, '_translated' if translated else '')
    df = pd.read_csv(url)
    return df
