from __future__ import annotations

from typing import Any, Callable, Literal

from ._db.omnipath import omnipath_universe

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


class Universe:


    def __init__(
            self,
            resource: Literal['omnipath'] | pd.DataFrame,
            **param
        ):
        """
        Load and preprocess a generic network from databases or files.
        """

        self._resource = resource
        self._param = param


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
