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
