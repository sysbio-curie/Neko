import pytest
import pandas as pd
from neko.inputs.signor import download_signor_database
from neko.inputs._db.signor import signor as process_signor_tsv
from neko.inputs._db.omnipath import omnipath_universe


def test_download_signor_database(tmp_path):
    # Download and save the SIGNOR dataset
    save_path = tmp_path / "SIGNOR_Human.tsv"
    download_signor_database(str(save_path))
    assert save_path.exists()
    # Load as DataFrame
    df = pd.read_csv(save_path, sep='\t')
    assert not df.empty
    # Process with the signor loader
    processed = process_signor_tsv(str(save_path))
    assert not processed.empty
    assert {'source', 'target'}.issubset(processed.columns)


def test_download_signor_database_as_df():
    # Download as DataFrame
    df = download_signor_database()
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert 'IDA' in df.columns and 'IDB' in df.columns


def test_omnipath_universe():
    # This test assumes omnipath is installed and available
    df = omnipath_universe()
    assert isinstance(df, pd.DataFrame)
    assert {'source', 'target'}.issubset(df.columns)
