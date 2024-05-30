from __future__ import annotations

from typing import Any, Callable, Iterable, Literal
import pickle
import logging
import pandas as pd

from pypath_common import misc as _common

from ._db.omnipath import omnipath_universe
from ..core import _networkbase as _nbase
from ._db import _misc

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

MANDATORY_COLUMNS = [
    'source',
    'target',
    'is_directed',
    'is_inhibition',
    'is_stimulation',
    'form_complex',
]

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
            resources: Literal['omnipath'] | pd.DataFrame = None,
            **param
        ):
        """
        Load and preprocess a generic network from databases or files.
        """

        self._resources = {}
        self.interactions = None
        self.add_resources(resources, **param)
        self.build()


    def add_resources(self, resources, **param) -> None:

        if isinstance(resources, str):

            if resources.endswith('.pickle'):

                with open(resources, 'rb') as fin:

                    resources = pickle.load(fin)

            elif resources.endswith('.tsv'):

                resources = pd.read_csv(resources, sep='\t')

        name = param.get('name', '_default')
        columns = param.get('columns', {})

        if isinstance(resources, pd.DataFrame):

            self._resources[name] = self._check_columns(resources, columns)

        elif isinstance(resources, str) and resources in _METHODS:

            name = resources

            self._resources[name] = self._check_columns(
                _METHODS[resources](**param),
                columns,
            )

        elif isinstance(resources, dict):

            resources = {
                k: self._check_columns(v, columns)
                for k, v in resources.items()
            }
            self._resources.update(resources)

        elif isinstance(resources, Iterable):

            for r in resources:

                self.add_resources(r, **param)


    @staticmethod
    def _check_columns(df: pd.DataFrame, columns: dict) -> pd.DataFrame:
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

        return df


    @staticmethod
    def merge(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
        """
        This function concatenates the provided df with the existing one in the resources object,
        aligning columns and filling in missing data with NaN.

        Parameters:
            df1 (pd.DataFrame): The DataFrame to be added.
            df2 (pd.DataFrame): The DataFrame to be added.

        Raises:
            ValueError: If the 'df' parameter is not a pandas DataFrame.

        Returns:
            None
        """

        # Align columns of both dataframes, filling missing columns with NaN
        all_columns = set(df1.columns).union(set(df2.columns))
        df1 = df1.reindex(columns=all_columns, fill_value=None)
        df2 = df2.reindex(columns=all_columns, fill_value=None)
        df1 = pd.concat([df1, df2])

        return df1.copy()


    def build(self, resources: Iterable[str] | None = None) -> None:

        resources = _common.to_list(resources) or self._resources.keys()

        self.interactions = None

        for res in resources:
            df = self._resources[res]
            self.interactions = (
                df.copy()
                if self.interactions is None else
                self.merge(self.interactions, df)
            )

        self.interactions.reset_index(drop=True, inplace=True)


    @property
    def resource(self):

        return self._resource if isinstance(self._resource, str) else 'user'


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

        return f'Universe; resources: {", ".join(self._resources.keys())}; size: {len(self)}'


    def __len__(self) -> int:

        return 0 if self.interactions is None else len(self.interactions)


    def check(self) -> bool:
        """
        The network is loaded and contains the mandatory variables.
        """

        return (
            hasattr(self, '_network') and
            not _REQUIRED_COLS - set(self._network.columns)
        )

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
