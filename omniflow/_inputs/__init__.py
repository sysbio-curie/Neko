from typing import TYPE_CHECKING

if TYPE_CHECKING:

    import pandas as pd

from ._db import omnipath_universe

"""
Access to generic networks from databases, files and standard formats.
"""


def network_universe(**kwargs) -> pd.DataFrame:
    """
    Generic networks from databases, files and standard formats.

    Args:
        kwargs:
            Passed to the source specific method. See the specific methods in
            this module for details.

    Note: currently OmniPath PPI is the single available option and serves
    as a placeholder. Later we will dispatch all inputs through this API.
    """

    return omnipath_universe(**kwargs)
