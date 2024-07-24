import pandas as pd


def signor(path: str) -> pd.DataFrame:
    """
    SIGNOR database from TSV.

    Processes SIGNOR interactions from a local TSV file.

    Parameters:
        path (str):
            The path to the SIGNOR TSV.

    Returns:
        None
    """

    # Read the SIGNOR database file into a pandas DataFrame
    df = pd.read_table(path)

    def contains(field, substr):

        return substr in field

    # Filter out the rows where EFFECT is "form complex" or "unknown"
    df = df[~df['EFFECT'] == "unknown"]

    # Transform the original dataframe into the desired format
    df = pd.DataFrame({
        'source': df['IDA'],
        'target': df['IDB'],
        'is_directed': df['DIRECT'],
        'is_stimulation': df['EFFECT'].apply(contains, 'up-regulates'),
        'is_inhibition': df['EFFECT'].apply(contains, 'down-regulates'),
        'form_complex': df['EFFECT'].apply(contains, 'complex'),
        'consensus_direction': False,  # Assuming no data provided, set all to False
        'consensus_stimulation': False,  # Assuming no data provided, set all to False
        'consensus_inhibition': False,  # Assuming no data provided, set all to False
        'curation_effort': df['ANNOTATOR'],
        'references': df['PMID'],
        'sources': df['SIGNOR_ID']
    })

    return df
