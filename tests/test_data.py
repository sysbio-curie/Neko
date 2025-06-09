import pandas as pd
import pytest
from neko.data import phosphosite_kinase_substrate, phosphosite_regulatory_sites

def test_phosphosite_kinase_substrate():
    df = phosphosite_kinase_substrate()
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert 'KINASE' in df.columns or 'Kinase' in df.columns

def test_phosphosite_regulatory_sites():
    df = phosphosite_regulatory_sites()
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert 'GENE' in df.columns or 'GENE' in df.columns

