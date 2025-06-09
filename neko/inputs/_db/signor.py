import pandas as pd
import requests
import pandas as pd
from io import BytesIO


def download_signor_database(save_path: str = None) -> pd.DataFrame:
    """
    Download the SIGNOR human dataset from the official API.
    If save_path is provided, saves the file as TSV. Otherwise, returns a DataFrame.
    """
    url = "https://signor.uniroma2.it/API/getHumanData.php"
    try:
        r = requests.get(url)
        r.raise_for_status()
        if save_path:
            with open(save_path, 'wb') as f:
                f.write(r.content)
            return None
        else:
            df = pd.read_csv(BytesIO(r.content), sep='\t')
            return df
    except requests.RequestException as e:
        raise RuntimeError(f"Error downloading SIGNOR database: {str(e)}")


def signor(path: str = None) -> pd.DataFrame:
    """
    SIGNOR database from TSV or download.

    Processes SIGNOR interactions from a local TSV file or downloads if no path is provided.

    Parameters:
        path (str, optional):
            The path to the SIGNOR TSV. If None, downloads the database.

    Returns:
        pd.DataFrame: Processed SIGNOR interactions.
    """
    if path is None:
        df = download_signor_database()
    else:
        df = pd.read_table(path)

    def contains(field, substr):
        return substr in field

    # Filter out the rows where EFFECT is "form complex" or "unknown"
    df = df[df['EFFECT'] != "unknown"]

    # Transform the original dataframe into the desired format
    df = pd.DataFrame({
        'source': df['IDA'],
        'target': df['IDB'],
        'is_directed': df['DIRECT'],
        'is_stimulation': df['EFFECT'].apply(lambda x: contains(x, 'up-regulates')),
        'is_inhibition': df['EFFECT'].apply(lambda x: contains(x, 'down-regulates')),
        'form_complex': df['EFFECT'].apply(lambda x: contains(x, 'complex')),
        'consensus_direction': False,  # Assuming no data provided, set all to False
        'consensus_stimulation': False,  # Assuming no data provided, set all to False
        'consensus_inhibition': False,  # Assuming no data provided, set all to False
        'curation_effort': df['ANNOTATOR'],
        'references': df['PMID'],
        'sources': df['SIGNOR_ID'],
    })

    # Add the transformed DataFrame to the existing database
    df = _group_by_source_target(df)

    return df


def _group_by_source_target(df_ungrouped: pd.DataFrame) -> pd.DataFrame:

    grouped = df_ungrouped.groupby(['source', 'target'])

    def join_strings(series):
        return '; '.join(set(series.astype(str)))

    def majority_true(series):
        # Count True values and compare to half the length of the series
        return series.sum() >= (len(series) / 2)

    def process_interaction(df):
        new_df = pd.DataFrame({
            'source': df['source'],
            'target': df['target'],
            'is_directed': df['is_directed'],
            'is_stimulation': df['is_stimulation'],
            'is_inhibition': df['is_inhibition'],
            'form_complex': df['form_complex'],
            'consensus_direction': True,
            'consensus_stimulation': (
                df['is_stimulation'].sum() >
                df['is_inhibition'].sum()
            ),
            'consensus_inhibition':  (
                df['is_inhibition'].sum() >
                df['is_stimulation'].sum()
            ),
            'curation_effort': df['curation_effort'],
            'references': df['references'],
            'sources': df['sources']
        })
        return new_df

    new_groups = []

    for (source, target), group_data in grouped:
        if len(group_data) > 1:
            aggregate_functions = {
                "is_stimulation": majority_true,
                "is_inhibition": majority_true,
                "form_complex": majority_true,
                "curation_effort": join_strings,
                "references": join_strings,
                "sources": join_strings
            }
            new_df = group_data.groupby(
                ["source", "target", "is_directed", "consensus_direction", "consensus_stimulation",
                 "consensus_inhibition"]).aggregate(aggregate_functions).reset_index()
            new_groups.append(process_interaction(new_df))
        else:
            new_groups.append(group_data)

    # Concatenate all groups back into a single DataFrame
    df_unique = pd.concat(new_groups, ignore_index=True)

    return df_unique

