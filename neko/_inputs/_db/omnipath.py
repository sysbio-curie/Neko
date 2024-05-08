import omnipath as op
import pandas as pd
import logging
import pypath
from pypath.utils import mapping

"""
Access to network databases.
"""


def omnipath_universe(**kwargs):
    """
    Access generic networks from OmniPath.
    """

    return op.interactions.PostTranslationalInteractions.get(**kwargs)


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
        self.interactions = op.interactions.PostTranslationalInteractions.get()
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

    def add_database(
            self,
            database: pd.DataFrame,
            column_mapping: dict = None,
            axis=0,
            ignore_index=True,
            reset_index=False,
        ):
        """
        This function concatenates the provided database with the existing one in the resources object,
        aligning columns and filling in missing data with NaN.

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

            # If only one of is_inhibition or is_activation is specified in the column_mapping, assign the values
            # from the specified column to both is_inhibition and is_activation, and then invert the values for one
            # of them
            if 'is_inhibition' in column_mapping and 'is_activation' not in column_mapping:
                if 'is_inhibition' in database.columns and database['is_inhibition'].nunique() > 1:
                    database['is_inhibition'] = database['is_inhibition'].replace(
                        {"True": True, "False": False, "true": True, "false": False, 1: True, -1: False}).astype(bool)
                    database['is_activation'] = ~database['is_inhibition']
            elif 'is_activation' in column_mapping and 'is_inhibition' not in column_mapping:
                if 'is_activation' in database.columns and database['is_activation'].nunique() > 1:
                    database['is_activation'] = database['is_activation'].replace(
                        {"True": True, "False": False, "true": True, "false": False, 1: True, -1: False}).astype(bool)
                    database['is_inhibition'] = ~database['is_activation']

        # Convert is_inhibition and is_activation columns to boolean values
        for column in ['is_inhibition', 'is_activation']:
            if column in database.columns:
                database[column] = database[column].replace({1: True, -1: False}).astype(bool)

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

    @property
    def nodes(self) -> set[str]:

        return (
            set(self.interactions['source']) |
            set(self.interactions['target'])
        )


    def __contains__(self, other) -> bool:

        return other in self.nodes


    def __and__(self, other: set) -> set:

        return self.nodes & other
