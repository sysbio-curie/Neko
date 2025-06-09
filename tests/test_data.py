import pandas as pd
import pytest
from neko.data import phosphosite_kinase_substrate, phosphosite_regulatory_sites, huri

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

def test_huri():
    df = huri()
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert 'GeneA' in df.columns or 'GeneB' in df.columns or 'A' in df.columns
