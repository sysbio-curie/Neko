import pandas as pd
from typing import Literal
def psp(kinase_substrate: str,
        regulatory_sites: str,
        organism: Literal["human", "mouse", "rat"],
        expand: bool = False
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
            expand (bool, optional): If True, expands each line to include a new interaction from phosphosite to
            respective protein. Default is False.

        Returns:
            ks (pd.DataFrame): The processed interactions dataframe.
    """

    # Load the kinase interaction and phosphorylation effect files
    ks = pd.read_csv(kinase_substrate, sep="\t")
    rs = pd.read_csv(regulatory_sites, sep="\t")

    # Filter the kinase interaction dataframe to include only interactions from the specified organism
    ks = ks.loc[
        (ks['KIN_ORGANISM'] == organism) &
        (ks['SUB_ORGANISM'] == organism)
    ]

    # Filter the phosphorylation effect dataframe to include only interactions from the specified organism
    rs = rs.loc[rs['ORGANISM'] == organism]

    # Concatenate the SUB_GENE and SUB_MOD_RSD columns to create the target column
    ks['target'] = ks['SUB_GENE'] + '_' + ks['SUB_MOD_RSD']

    # Create the ks dataframe with the source and target columns
    ks = pd.DataFrame({
        'source': ks['GENE'],
        'target': ks['target'],
    })

    # Keep only MOD_RSD entries with "-p" suffix and remove it
    rs = rs[rs['MOD_RSD'] == rs['MOD_RSD'].endswith('-p')]
    rs['MOD_RSD'] = rs['MOD_RSD'].apply(lambda x: x[:-2])
    # Concatenate the GENE and MOD_RSD columns to create the Prot_site column
    rs['target'] = rs['GENE'] + '_' + rs['MOD_RSD']

    rs['is_stimulation'] = rs['ON_FUNCTION'].str.contains('activity, induced')
    rs['is_inhibition'] = rs['ON_FUNCTION'].str.contains('activity, inhibited')

    # Update the ks dataframe based on the phosphorylation effect information
    ks = ks.merge(rs[['target', 'is_inhibition', 'is_stimulation']], how='left')
    ks.fillna('', inplace=True)

    # If expand is True, expand each line to include a new interaction from phosphosite to respective protein
    if expand:
        ptm_protein = rs.copy()
        ptm_protein['source'] = ptm_protein['target']
        ptm_protein['target'] = ptm_protein['target'].apply(lambda x: x.split('_')[0])
        ks = pd.concat([ks, pd.DataFrame(ptm_protein)], ignore_index=True)

    # Filter the ks dataframe based on the is_stimulation and is_inhibition columns
    ks = ks[ks['is_stimulation'] | ks['is_inhibition']]

    # If expand is False, remove duplicates from the dataframe
    if not expand:
        ks = ks.drop_duplicates(subset=['source', 'target'])

    return ks
