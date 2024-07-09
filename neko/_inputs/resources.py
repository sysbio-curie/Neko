import omnipath as op
import pandas as pd
import logging
import pypath
from pypath.utils import mapping


class Resources():
    """
    This class is used to store and manage databases for mining interactions.
    The user can select different database formats from omnipath, pypath or load one of their own.

    Attributes:
    interactions (pd.DataFrame): The main DataFrame that stores the interactions.
    filtered_interactions (pd.DataFrame): A filtered version of the interactions DataFrame. Currently not implemented.
    required_columns (list): A list of column names that are required in the interactions DataFrame.
    """

    def __init__(self):
        """
        The constructor for the Resources class. Initializes the interactions and filtered_interactions attributes as None.
        Sets up logging configuration and defines the required columns for the interactions DataFrame.
        """

        # The main DataFrame that stores the interactions
        self.interactions = None

        # A filtered version of the interactions DataFrame. Currently not implemented
        self.filtered_interactions = None

        # Set up logging configuration
        logging.basicConfig(level=logging.INFO)

        # Define the required columns for the interactions DataFrame
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

    def add_database(self, database: pd.DataFrame, column_mapping: dict = None, axis=0, ignore_index=True,
                     reset_index=False):
        """
        This function concatenates the provided database with the existing one in the resources object,
        aligning columns and filling in missing data with NaN. It also ensures that all required columns
        are present, creating them with default values if necessary.

        Parameters:
        database (pd.DataFrame): The DataFrame to be added.
        column_mapping (dict, optional): A dictionary mapping column names in the incoming database to column names in the resources object.
        axis (int, optional): 0 to concatenate row-wise, 1 to concatenate column-wise. Default is 0.
        ignore_index (bool, optional): If True, the resulting axis will be labeled 0, 1, …, n - 1. Default is True.
        reset_index (bool, optional): If True, reset the index of the concatenated DataFrame. Default is False.

        Raises:
        ValueError: If the 'database' parameter is not a pandas DataFrame.

        Returns:
        None
        """

        # Check if the provided database is a pandas DataFrame
        if not isinstance(database, pd.DataFrame):
            raise ValueError("The 'database' parameter must be a pandas DataFrame.")

        # If column_mapping is provided, rename the columns of the incoming database
        if column_mapping is not None:
            database = database.rename(columns=column_mapping)

        # Ensure all required columns are present, create if missing with default values
        for column in self.required_columns:
            if column not in database.columns:
                if column in ['is_directed', 'is_stimulation', 'is_inhibition', 'form_complex',
                              'consensus_direction', 'consensus_stimulation', 'consensus_inhibition']:
                    database[column] = False
                elif column in ['curation_effort']:
                    database[column] = None
                elif column in ['references', 'sources']:
                    database[column] = ''
                else:
                    database[column] = None
                logging.info(f"Created missing column '{column}' with default values.")

        # Handle is_inhibition and is_activation columns
        if 'is_inhibition' in database.columns and 'is_activation' not in database.columns:
            database['is_inhibition'] = database['is_inhibition'].replace(
                {"True": True, "False": False, "true": True, "false": False, 1: True, -1: False}).astype(bool)
            database['is_activation'] = ~database['is_inhibition']
        elif 'is_activation' in database.columns and 'is_inhibition' not in database.columns:
            database['is_activation'] = database['is_activation'].replace(
                {"True": True, "False": False, "true": True, "false": False, 1: True, -1: False}).astype(bool)
            database['is_inhibition'] = ~database['is_activation']

        # Convert boolean columns to proper boolean values
        bool_columns = ['is_directed', 'is_stimulation', 'is_inhibition', 'form_complex',
                        'consensus_direction', 'consensus_stimulation', 'consensus_inhibition']
        for column in bool_columns:
            if column in database.columns:
                database[column] = database[column].replace({1: True, 0: False, -1: False}).astype(bool)

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

        logging.info("Database added successfully.")
        return

    def import_signor_tsv(self, signor_file):
        """
        This function imports a SIGNOR database file in TSV format, processes it, and adds it to the existing database.

        Parameters:
        signor_file (str): The path to the SIGNOR database file in TSV format.

        Returns:
        None
        """

        # Read the SIGNOR database file into a pandas DataFrame
        df_signor = pd.read_table(signor_file)

        # Define a function to determine if the effect is stimulation
        def is_stimulation(effect):
            """
            This function checks if the effect is stimulation.

            Parameters:
            effect (str): The effect to check.

            Returns:
            bool: True if the effect is stimulation, False otherwise.
            """
            return 'up-regulates' in effect

        # Define a function to determine if the effect is inhibition
        def is_inhibition(effect):
            """
            This function checks if the effect is inhibition.

            Parameters:
            effect (str): The effect to check.

            Returns:
            bool: True if the effect is inhibition, False otherwise.
            """
            return 'down-regulates' in effect

        # Define a function to determine if the effect is form complex
        def form_complex(effect):
            """
            This function checks if the effect is form complex.

            Parameters:
            effect (str): The effect to check.

            Returns:
            bool: True if the effect is form complex, False otherwise.
            """
            return 'complex' in effect

        def is_directed(directed):
            """
            This function checks if the interaction is directed.

            Parameters:
            directed (str): The value to check.

            Returns:
            bool: True if the interaction is directed, False otherwise.
            """
            return True if directed == 'YES' else False

        # Filter out the rows where EFFECT is "form complex" or "unknown"
        filtered_df = df_signor[~df_signor['EFFECT'].isin(["unknown"])]

        # Transform the original dataframe into the desired format
        transformed_df = pd.DataFrame({
            'source': filtered_df['IDA'],
            'target': filtered_df['IDB'],
            'is_directed': filtered_df['DIRECT'].apply(is_directed),
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

        # Add the transformed DataFrame to the existing database
        grouped_database = self.group_by_source_target(transformed_df)
        self.add_database(grouped_database)

        return

    def group_by_source_target(self, df_ungrouped):
        grouped = df_ungrouped.groupby(['source', 'target'])

        def join_strings(series):
            return '; '.join(set(series.astype(str)))

        def majority_true(series):
            # Count True values and compare to half the length of the series
            return series.sum() >= (len(series) / 2)

        def process_consensus(df):
            new_df = pd.DataFrame({
                'source': df['source'],
                'target': df['target'],
                'is_directed': df['is_directed'],
                'is_stimulation': df['is_stimulation'],
                'is_inhibition': df['is_inhibition'],
                'form_complex': df['form_complex'],
                'consensus_direction': True,
                'consensus_stimulation': df['is_stimulation'],
                'consensus_inhibition': df['is_inhibition'],
                'curation_effort': df['curation_effort'],
                'references': df['references'],
                'sources': df['sources']
            })
            return new_df

        def process_interaction(df):
            new_df = pd.DataFrame({
                'source': df['source'],
                'target': df['target'],
                'is_directed': df['is_directed'],
                'is_stimulation': df['is_stimulation'],
                'is_inhibition': df['is_inhibition'],
                'form_complex': df['form_complex'],
                'consensus_direction': True,
                'consensus_stimulation': False,
                'consensus_inhibition': False,
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
                if len(group_data[group_data["is_stimulation"] == True]) == len(
                    group_data[group_data["is_inhibition"] == True]):
                    new_groups.append(process_interaction(new_df))
                else:
                    new_groups.append(process_consensus(new_df))
            else:
                new_groups.append(group_data)

        # Concatenate all groups back into a single DataFrame
        df_unique = pd.concat(new_groups, ignore_index=True)

        return df_unique

    def process_psp_interactions(self, kinase_int_file, phospho_effect_file, organism, expand=False):
        """
        This function loads files from PhoshpositePlus (PSP), parses the files to create sign interactions based on the effect of phosphorylation on protein activities
        and creates an interaction dataframe based on the Omnipath interaction format.

        Parameters:
        kinase_int_file (str): The path to the kinase interaction file.
        phospho_effect_file (str): The path to the phosphorylation effect file.
        organism (str): The organism to filter interactions for.
        expand (bool, optional): If True, expands each line to include a new interaction from phosphosite to respective protein. Default is False.

        Returns:
        psp_interactions (pd.DataFrame): The processed interactions dataframe.
        """

        # Load the kinase interaction and phosphorylation effect files
        kinase_int = pd.read_csv(kinase_int_file, sep="\t")
        phospho_effect = pd.read_csv(phospho_effect_file, sep="\t")

        # Filter the kinase interaction dataframe to include only interactions from the specified organism
        kinase_int_filtered = kinase_int.loc[
            (kinase_int['KIN_ORGANISM'] == organism) & (kinase_int['SUB_ORGANISM'] == organism)]

        # Concatenate the SUB_GENE and SUB_MOD_RSD columns to create the target column
        kinase_int_filtered['target'] = kinase_int_filtered['SUB_GENE'] + '_' + kinase_int_filtered['SUB_MOD_RSD']

        # Create the psp_interactions dataframe with the source and target columns
        psp_interactions = pd.DataFrame({
            'source': kinase_int_filtered['GENE'],
            'target': kinase_int_filtered['target'],
            'is_directed': True,
            'consensus_direction': False,  # Assuming no data provided, set all to False
            'consensus_stimulation': False,  # Assuming no data provided, set all to False
            'consensus_inhibition': False,  # Assuming no data provided, set all to False
        })

        # Filter the phosphorylation effect dataframe to include only interactions from the specified organism
        phospho_effect_filtered = phospho_effect.loc[phospho_effect['ORGANISM'] == organism]

        # Keep only MOD_RSD entries with "-p" suffix and remove it
        phospho_effect_filtered['MOD_RSD'] = phospho_effect_filtered['MOD_RSD'].apply(
            lambda x: x[:-2] if x.endswith('-p') else x)

        # Concatenate the GENE and MOD_RSD columns to create the Prot_site column
        phospho_effect_filtered['Prot_site'] = phospho_effect_filtered['GENE'] + '_' + phospho_effect_filtered[
            'MOD_RSD']

        # Initialize the is_stimulation and is_inhibition columns with 0
        psp_interactions['is_stimulation'] = 0
        psp_interactions['is_inhibition'] = 0

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

        # Update the psp_interactions dataframe based on the phosphorylation effect information
        psp_interactions = psp_interactions.merge(phospho_effect_filtered[['Prot_site', 'ON_FUNCTION']],
                                                  left_on='target', right_on='Prot_site', how='left')
        psp_interactions.fillna('', inplace=True)
        psp_interactions = psp_interactions.apply(update_phospho_effect, axis=1)
        psp_interactions.drop(columns=['Prot_site', 'ON_FUNCTION'], inplace=True)

        # If expand is True, expand each line to include a new interaction from phosphosite to respective protein
        if expand:
            expanded_rows = []
            for index, row in psp_interactions.iterrows():
                source, target = row['source'], row['target']
                target_protein = target.split('_')[0]
                expanded_rows.append({'source': target, 'target': target_protein, 'is_directed': True,
                                      'is_stimulation': row['is_stimulation'], 'is_inhibition': row['is_inhibition'],
                                      'consensus_direction': False, 'consensus_stimulation': False,
                                      'consensus_inhibition': False, 'curation_effort': False, 'references': False,
                                      'sources': False})
            psp_interactions = pd.concat([psp_interactions, pd.DataFrame(expanded_rows)], ignore_index=True)

        # Filter the psp_interactions dataframe based on the is_stimulation and is_inhibition columns
        psp_interactions = psp_interactions[
            (psp_interactions['is_stimulation'] != 0) | (psp_interactions['is_inhibition'] != 0)]

        # If expand is False, remove duplicates from the dataframe
        if not expand:
            psp_interactions = psp_interactions.drop_duplicates(subset=['source', 'target'])

        # Add the processed interactions to the existing database
        self.add_database(psp_interactions)

        return psp_interactions
