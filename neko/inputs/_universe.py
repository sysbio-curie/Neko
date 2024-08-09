from __future__ import annotations

from typing import Any, Callable, Iterable, Literal
import pickle
import logging
import pandas as pd
import os
from pypath_common import misc as _common

from ._db.omnipath import omnipath_universe
from ._db import psp as _psp
from ._db import _misc
from ._db import huri as _huri
from ._db import signor as _signor

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

MANDATORY_BOOL_COLS = [
    'is_directed',
    'is_stimulation',
    'is_inhibition',
    'form_complex',
]

MANDATORY_COLUMNS = [
    'source',
    'target',
]
MANDATORY_COLUMNS.extend(MANDATORY_BOOL_COLS)


def network_universe(
        resource: Literal['omnipath'] | pd.DataFrame = 'omnipath',
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


def omnipath(**kwargs) -> Universe:

    return network_universe('omnipath', **kwargs)


def signor(path: str | None = None, **kwargs) -> Universe:

    if path and os.path.exists(path):

        return Universe(_signor.signor(path))

    else:

        raise NotImplementedError(
            'Direct retrieval of SIGNOR is not implemented. '
            'Please download manually the TSV file and provide the path.'
        )


def phosphosite(
        organism: Literal["human", "mouse", "rat"] = 'human',
        kinase_substrate: str | None = None,
        regulatory_sites: str | None = None,
        **kwargs
    ) -> Universe:

    df = _psp.psp(organism, kinase_substrate, regulatory_sites)

    return Universe(df, name = 'phosphosite')


def huri(dataset: str = 'HI-union') -> Universe:

    df = _huri.huri(dataset)

    return Universe(df, directed = False, name = 'huri')


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
        self._directed = {}
        self.interactions = None
        self.add_resources(resources, **param)
        self.build()


    def add_resources(self, resources, directed: bool = True, **param) -> None:

        if isinstance(resources, str):

            if resources.endswith('.pickle'):

                with open(resources, 'rb') as fin:

                    resources = pickle.load(fin)

            elif resources.endswith('.tsv'):

                resources = pd.read_csv(resources, sep='\t')

        name = param.get('name', '_default')
        columns = param.get('columns', {})

        if isinstance(resources, str) and resources in _METHODS:

            name = name or resources
            resources = _METHODS[resources](**param)

        if isinstance(resources, pd.DataFrame):

            self._resources[name] = self._check_columns(
                resources,
                columns,
                directed = directed,
            )

        elif isinstance(resources, dict):

            resources = {
                k: self._check_columns(v, columns)
                for k, v in resources.items()
            }

            if isinstance(directed, bool):

                directed = {k: directed for k in resources.keys()}

            for name, resource in resources.items():

                self.add_resources(
                    resource,
                    directed = directed[name],
                    name = name,
                    **param,
                )

        elif isinstance(resources, Iterable):

            names = [f'{name}_{i}' for i in range(len(resources))]

            if isinstance(directed, bool):

                directed = [directed] * len(resources)

            for r, d, n in zip(resources, directed, names):

                self.add_resources(r, directed = d, name = n, **param)


    @staticmethod
    def _check_columns(
            df: pd.DataFrame,
            columns: dict,
            directed: bool = True,
        ) -> pd.DataFrame:

        # If columns is provided, rename the columns of the incoming df
        if columns:
            df = df.rename(columns=columns)

        if 'effect' in df.columns:

            df = _misc.split_effect(df)

        for col in MANDATORY_BOOL_COLS:

            df = _misc.bool_col(df, col)

        # Check if the df contains the required columns
        missing_columns = set(MANDATORY_COLUMNS) - set(df.columns)

        if missing_columns:

            logging.warning("The incoming df is missing some required columns: %s", missing_columns)
            logging.warning("This might lead to issues in running the package.")

        if not directed:

            df = _misc.undirected_to_mutual(df)

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

        if not resources:

            return

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
