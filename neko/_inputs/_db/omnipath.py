import omnipath as op
import pandas as pd
import logging
import pypath
from pypath.utils import mapping

import _misc as _misc

"""
Access to network databases.
"""


def omnipath_universe(**kwargs):
    """
    Access generic networks from OmniPath.
    """

    return op.interactions.PostTranslationalInteractions.get(**kwargs)


MANDATORY_COLUMNS = [
    'source',
    'target',
    'is_directed',
    'is_inhibition',
    'is_stimulation',
    'form_complex',
]


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
            df: pd.DataFrame,
            columns: dict = None,
        ):
        """
        This function concatenates the provided df with the existing one in the resources object,
        aligning columns and filling in missing data with NaN.

        Parameters:
            df (pd.DataFrame): The DataFrame to be added.
            columns (dict, optional):
                A dictionary of column name mappings to be applied on the
                provided data frame. Mandatory columns: source, target, is_directed,
                is_inhibition, is_stimulation, form_complex.

        Raises:
            ValueError: If the 'df' parameter is not a pandas DataFrame.

        Returns:
            None
        """

        # If columns is provided, rename the columns of the incoming df
        if columns:
            df = df.rename(columns=columns)

        if 'effect' in df.columns:

            df = _misc.split_effect(df)


        df = _misc.bool_col(df, 'is_stimulation')
        df = _misc.bool_col(df, 'is_inhibition')

        # Check if the df contains the required columns
        missing_columns = set(MANDATORY_COLUMNS) - set(df.columns)

        if missing_columns:

            logging.warning("The incoming df is missing some required columns: %s", missing_columns)
            logging.warning("This might lead to issues in running the package.")

        if self.interactions is not None:

            # Align columns of both dataframes, filling missing columns with NaN
            all_columns = set(self.interactions.columns).union(set(df.columns))
            self.interactions = self.interactions.reindex(columns=all_columns, fill_value=None)
            df = df.reindex(columns=all_columns, fill_value=None)
            self.interactions = pd.concat([self.interactions, df])

        elif self.interactions is None:
            # If self.interactions is None, initialize it with the incoming df
            self.interactions = df

        self.interactions.reset_index(drop=True, inplace=True)


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
