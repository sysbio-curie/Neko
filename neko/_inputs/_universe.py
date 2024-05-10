from __future__ import annotations

from typing import Any, Callable, Literal
import pickle
import pandas as pd

from ._db.omnipath import omnipath_universe
from ..core import _networkbase as _nbase

"""
Access to generic networks from databases, files and standard formats.
"""

_METHODS = {
    'omnipath': omnipath_universe,
}
_REQUIRED_COLS = {
    'source',
    'target',
}


def network_universe(
        resource: Literal['omnipath'] | pd.DataFrame,
        **kwargs
    ) -> Universe:
    """
    Generic networks from databases, files and standard formats.

    Args:
        resource:
            Name of the resource or a ready data frame to bypass the built-in
            loading method.
        kwargs:
            Passed to the source specific method. See the specific methods in
            this module for details.

    Note: currently OmniPath PPI is the single available option and serves
    as a placeholder. Later we will dispatch all inputs through this API.
    """

    return Universe(resource, **kwargs)


def omnipath(**kwargs):

    return network_universe('omnipath', **kwargs)


class Universe:


    def __init__(
            self,
            resources: Literal['omnipath'] | pd.DataFrame,
            **param
        ):
        """
        Load and preprocess a generic network from databases or files.
        """

        self.add_resources(resources, **param)


    def add_resources(self, resources, **param) -> None:

        if isinstance(resources, str):

            if resources.endswith('.pickle'):

                with open(resources, 'rb') as fin:

                    resources = pickle.load(fin)

            elif resources.endswith('.tsv'):

                resources = pd.read_csv(resources, sep='\t')

        if isinstance(resources, pd.DataFrame):

            name = param.get('name', '_default')
            self._resources[name] = resources

        elif isinstance(resources, str) and resources in _METHODS:

            self._resources[name] = _METHODS[resources](**param)

        elif isinstance(resources, dict):

            self._resources.update(resources)

        elif isinstance(resources, Iterable):

            for r in resources:

                self.add_resources(r, **param)





    @staticmethod
    def merge(
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
    def resource(self):

        return self._resource if isinstence(self._resource, str) else 'user'


    @property
    def network(self) -> pd.DataFrame:
        """
        The network as it's been read from the original source.
        """

        if not hasattr(self, '_network'):

            self.load()

        return self._network


    @network.setter
    def network(self, value: Any):

       raise AttributeError('The attribute `Universe.network` is read-only.')


    def load(self) -> None:
        """
        Acquire the input data according to parameters.
        """

        self._network = self.method(**self.param)


    @property
    def method(self) -> Callable:
        """
        The method that loads the data.

        ``param`` are to be passed to this method.
        """

        return (
            lambda **kwargs: self._resource
                if isinstance(self._resource, pd.DataFrame) else
            _METHODS[self.resource]
        )


    def __repr__(self) -> str:

        param = (
            _common.dict_str(self._param)[42:].
            rsplit(', ', maxsplit = 1)[0]
        )
        param = f'; [{param}]' if param else ''

        return f'Universe from {self.resource}; size: {len(self)}{param}'


    def __len__(self) -> int:

        return len(getattr(self, '_network', ()))


    def check(self) -> bool:
        """
        The network is loaded and contains the mandatory variables.
        """

        return (
            hasattr(self, '_network') and
            not _REQUIRED_COLS - set(self._network.columns)
        )


MANDATORY_COLUMNS = [
    'source',
    'target',
    'is_directed',
    'is_inhibition',
    'is_stimulation',
    'form_complex',
]


class Resources(_nbase.NetworkBase):
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
