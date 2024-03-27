import omnipath as op
import pandas as pd
import logging
import pypath
from pypath.utils import mapping

class Resources():
    """
    This class stores the actual databases to mine interesting interactions. The user can select different
    database format from omnipath, pypath or loading one of its own.
    """

    def __init__(self):
        self.interactions = None
        # not sure the following will ever be implemented... maybe one day
        self.filtered_interactions = None
        # Configure logging as needed
        logging.basicConfig(level=logging.INFO)
        # Define the required columns
        self.required_columns = ['source', 'target', 'is_directed', 'is_stimulation', 'is_inhibition',
                                 'form_complex',
                                 'consensus_direction', 'consensus_stimulation', 'consensus_inhibition',
                                 'curation_effort', 'references', 'sources']

    def load_all_omnipath_interactions(self):
        """
        loads into the Resources object the omnipath dataframe "all_interactions"
        """
        self.interactions = op.interactions.AllInteractions.get()
        return

    def translate_dataframe_from_uniprot_to_genesymbol(self):
        """If the loaded dataframe has the source and target columns in uniprot, it will translate it to genesymbol"""

        return

    def load_omnipath_interactions(self):
        """
        loads into the Resources object the omnipath dataframe "Omnipath"
        """
        self.interactions = op.interactions.OmniPath.get()
        return

    def add_database(self, database: pd.DataFrame, axis=0, ignore_index=True, reset_index=False):
        """
        Concatenate the database that the user wants to add with the one in the resources object,
        aligning columns and filling in missing data with NaN.

        :param database: DataFrame to be added.
        :param axis: 0 to concatenate row-wise, 1 to concatenate column-wise. Default is 0.
        :param ignore_index: If True, the resulting axis will be labeled 0, 1, …, n - 1. Default is True.
        :param reset_index: If True, reset the index of the concatenated DataFrame. Default is False.
        """
        if not isinstance(database, pd.DataFrame):
            raise ValueError("The 'database' parameter must be a pandas DataFrame.")

        # Check if the database contains the required columns
        missing_columns = set(self.required_columns) - set(database.columns)
        if missing_columns:
            logging.warning("The incoming database is missing some required columns: %s", missing_columns)
            logging.warning("This might lead to issues in running the package.")

        # Align columns of both dataframes, filling missing columns with NaN
        if self.interactions is not None:
            all_columns = set(self.interactions.columns).union(set(database.columns))
            self.interactions = self.interactions.reindex(columns=all_columns, fill_value=None)
            database = database.reindex(columns=all_columns, fill_value=None)

        # If self.interactions is None, initialize it with the incoming database
        if self.interactions is None:
            self.interactions = database
        else:
            self.interactions = pd.concat([self.interactions, database], axis=axis, ignore_index=ignore_index)
            if reset_index:
                self.interactions.reset_index(drop=True, inplace=True)
        return

    def import_signor_tsv(self, signor_file):
        df_signor = pd.read_table(signor_file)
        # Substitute "up-regulates" with "stimulation" and "down-regulates" with "inhibition" in the EFFECT column

        # Function to determine if the effect is stimulation
        def is_stimulation(effect):
            return 'up-regulates' in effect

        # Function to determine if the effect is inhibition
        def is_inhibition(effect):
            return 'down-regulates' in effect

        # Function to determine if the effect is form complex
        def form_complex(effect):
            return 'complex' in effect

        # First, filter out the rows where EFFECT is "form complex" or "unknown"
        filtered_df = df_signor[~df_signor['EFFECT'].isin(["unknown"])]

        # Transform the original dataframe
        transformed_df = pd.DataFrame({
            'source': filtered_df['IDA'],
            'target': filtered_df['IDB'],
            'is_directed': filtered_df['DIRECT'],
            'is_stimulation': filtered_df['EFFECT'].apply(is_stimulation),
            'is_inhibition': filtered_df['EFFECT'].apply(is_inhibition),
            'form_complex': filtered_df['EFFECT'].apply(form_complex),
            'consensus_direction': False,  # Assuming no data provided, set all to False
            'consensus_stimulation': False,  # Assuming no data provided, set all to False
            'consensus_inhibition': False,  # Assuming no data provided, set all to False
            'curation_effort': filtered_df['ANNOTATOR'],
            'references': filtered_df['PMID'],
            'sources': filtered_df['SIGNOR_ID']
        })
        self.add_database(transformed_df)

        return

    def process_psp_interactions(self, kinase_int_file, phospho_effect_file, organism, expand=False):
        """
        Loads files from PhoshpositePlus (PSP), parses the files to create sign interactions based on the effect of phosphorylation on protein activities
        and creates an interaction dataframe based on the Omnipath interaction format.
        """
        kinase_int = pd.read_csv(kinase_int_file)
        phospho_effect = pd.read_csv(phospho_effect_file)

        # Filter kinase_int dataframe to include only interactions from specified organisms
        kinase_int_filtered = kinase_int.loc[(kinase_int['KIN_ORGANISM'] == organism) & (kinase_int['SUB_ORGANISM'] == organism)]

        # Concatenate SUB_GENE and SUB_MOD_RSD columns separated by "_"
        kinase_int_filtered['target'] = kinase_int_filtered['SUB_GENE'] + '_' + kinase_int_filtered['SUB_MOD_RSD']

        # Create psp_interactions dataframe with source and target columns
        psp_interactions = pd.DataFrame({
            'source': kinase_int_filtered['GENE'],
            'target': kinase_int_filtered['target'],
            'is_directed': True,
            'consensus_direction': False,  # Assuming no data provided, set all to False
            'consensus_stimulation': False,  # Assuming no data provided, set all to False
            'consensus_inhibition': False,  # Assuming no data provided, set all to False
        })

        # Filter phospho_effect dataframe to include only interactions from specified organisms
        phospho_effect_filtered = phospho_effect.loc[phospho_effect['ORGANISM'] == organism]

        # Keep only MOD_RSD entries with "-p" suffix and remove it
        phospho_effect_filtered['MOD_RSD'] = phospho_effect_filtered['MOD_RSD'].apply(lambda x: x[:-2] if x.endswith('-p') else x)

        # Concatenate GENE and MOD_RSD columns into a new column 'Prot_site'
        phospho_effect_filtered['Prot_site'] = phospho_effect_filtered['GENE'] + '_' + phospho_effect_filtered['MOD_RSD']

        # Initialize is_stimulation and is_inhibition columns with 0
        psp_interactions['is_stimulation'] = 0
        psp_interactions['is_inhibition'] = 0

        # Function to update is_stimulation and is_inhibition columns based on ON_FUNCTION values
        def update_phospho_effect(row):
            if isinstance(row['ON_FUNCTION'], str):  # Check if the value is a string
                activities = []
                # Split each value in ON_FUNCTION column
                for function in row['ON_FUNCTION'].split(';'):
                    # Check if the function contains the word "activity"
                    if 'activity' in function:
                        activities.append(function.strip())  # Add the function to the list of activities
                # Initialize stimulation and inhibition flags
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

        # Update psp_interactions dataframe based on phospho_effect information
        psp_interactions = psp_interactions.merge(phospho_effect_filtered[['Prot_site', 'ON_FUNCTION']], left_on='target', right_on='Prot_site', how='left')
        psp_interactions.fillna('', inplace=True)
        psp_interactions = psp_interactions.apply(update_phospho_effect, axis=1)
        psp_interactions.drop(columns=['Prot_site', 'ON_FUNCTION'], inplace=True)

        # If expand is True, expand each line to include a new interaction from phosphosite to respective protein
        if expand:
            expanded_rows = []
            for index, row in psp_interactions.iterrows():
                source, target = row['source'], row['target']
                target_protein = target.split('_')[0]
                expanded_rows.append({'source': target, 'target': target_protein, 'is_directed': True, 'is_stimulation': row['is_stimulation'], 'is_inhibition': row['is_inhibition'], 'consensus_direction': False, 'consensus_stimulation': False, 'consensus_inhibition': False, 'curation_effort': False, 'references': False, 'sources': False})
            psp_interactions = pd.concat([psp_interactions, pd.DataFrame(expanded_rows)], ignore_index=True)

            # psp_interactions dataframe based on phospho_effect information
        psp_interactions = psp_interactions[(psp_interactions['is_stimulation'] != 0) | (psp_interactions['is_inhibition'] != 0)]

        # If expand is False, remove duplicates from the DataFrame
        if not expand:
            psp_interactions = psp_interactions.drop_duplicates(subset=['source', 'target'])

        self.add_database(psp_interactions)
        return psp_interactions
