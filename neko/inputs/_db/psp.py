import pandas as pd
from typing import Literal

from ... import data as _data


def psp(
    organism: Literal["human", "mouse", "rat"],
    kinase_substrate: str | None = None,
    regulatory_sites: str | None = None,
) -> pd.DataFrame:
    """
    PhoshpositePlus (PSP) from local file.

    Parses the files to create sign interactions based on the
    effect of phosphorylation on protein activities and creates an interaction dataframe based on the Omnipath
    interaction format.

    Parameters:
        kinase_substrate (str): The path to the kinase interaction file.
        regulatory_sites (str): The path to the phosphorylation effect file.
        organism (str): The organism to filter interactions for.

    Returns:
        ks (pd.DataFrame): The processed interactions dataframe.
    """

    args = locals()

    def load(dset):

        nonlocal args
        value = args[dset]

        if isinstance(value, str):

            value = pd.read_table(value, sep='\t')

        elif value is None:

            value = getattr(_data, f'phosphosite_{dset}')()

        if not isinstance(value, pd.DataFrame):
            raise ValueError(f"Invalid type for '{dset}' ({type(value)}).")

        return value

    ks, rs = tuple(load(dset) for dset in ('kinase_substrate', 'regulatory_sites'))

    # Filter the kinase interaction dataframe to include only interactions from the specified organism
    ks = ks.loc[
        (ks['KIN_ORGANISM'] == organism) &
        (ks['SUB_ORGANISM'] == organism)
    ]

    # Concatenate the SUB_GENE and SUB_MOD_RSD columns to create the target column
    ks['target'] = ks['SUB_GENE'] + '_' + ks['SUB_MOD_RSD']

    # Create the ks dataframe with the source and target columns
    ks = pd.DataFrame({
        'source': ks['GENE'],
        'target': ks['target'],
        'is_directed': True,
        'consensus_direction': False,  # Assuming no data provided, set all to False
        'consensus_stimulation': False,  # Assuming no data provided, set all to False
        'consensus_inhibition': False,  # Assuming no data provided, set all to False
    })

    # Filter the phosphorylation effect dataframe to include only interactions from the specified organism
    rs = rs.loc[rs['ORGANISM'] == organism]

    # Keep only MOD_RSD entries with "-p" suffix and remove it
    rs = rs[rs['MOD_RSD'].str.endswith('-p')]
    rs['MOD_RSD'] = rs['MOD_RSD'].apply(lambda x: x[:-2])

    # Concatenate the GENE and MOD_RSD columns to create the Prot_site column
    rs['Prot_site'] = rs['GENE'] + '_' + rs[
        'MOD_RSD']

    # Initialize the is_stimulation and is_inhibition columns with 0
    ks['is_stimulation'] = 0
    ks['is_inhibition'] = 0

    # Function to update the is_stimulation and is_inhibition columns based on the ON_FUNCTION values
    def update_phospho_effect(row):
        """
        This function updates the is_stimulation and is_inhibition columns based on the ON_FUNCTION values.

        Parameters:
        row (pd.Series): The row to update.

        Returns:
        row (pd.Series): The updated row.
        """
        if isinstance(row['ON_FUNCTION'], str):  # Check if the value is a string
            activities = []
            # Split each value in the ON_FUNCTION column
            for function in row['ON_FUNCTION'].split(';'):
                # Check if the function contains the word "activity"
                if 'activity' in function:
                    activities.append(function.strip())  # Add the function to the list of activities
            # Initialize the stimulation and inhibition flags
            is_stimulation = 0
            is_inhibition = 0
            # Check each activity for "induced" or "inhibited"
            for activity in activities:
                if 'induced' in activity:
                    is_stimulation = 1
                elif 'inhibited' in activity:
                    is_inhibition = 1
            # Update the row with the flags
            row['is_stimulation'] = is_stimulation
            row['is_inhibition'] = is_inhibition

        return row

    # Update the ks dataframe based on the phosphorylation effect information
    ks = ks.merge(rs[['Prot_site', 'ON_FUNCTION']],
                  left_on='target', right_on='Prot_site', how='left')
    ks.fillna('', inplace=True)
    ks = ks.apply(update_phospho_effect, axis=1)
    ks.drop(columns=['Prot_site', 'ON_FUNCTION'], inplace=True)

    # If expand is True, expand each line to include a new interaction from phosphosite to respective protein
    expanded_rows = []

    for index, row in ks.iterrows():

        source, target = row['source'], row['target']
        target_protein = target.split('_')[0]
        expanded_rows.append({'source': target, 'target': target_protein, 'is_directed': True,
                              'is_stimulation': row['is_stimulation'], 'is_inhibition': row['is_inhibition'],
                              'consensus_direction': False, 'consensus_stimulation': False,
                              'consensus_inhibition': False, 'curation_effort': False, 'references': False,
                              'sources': False})

    ks = pd.concat([ks, pd.DataFrame(expanded_rows)], ignore_index=True)

    # Filter the ks dataframe based on the is_stimulation and is_inhibition columns
    ks = ks[
        (ks['is_stimulation'] != 0) | (ks['is_inhibition'] != 0)]

    ks = ks.drop_duplicates(subset=['source', 'target'])

    return ks
