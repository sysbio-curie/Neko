import omnipath as op
import pandas as pd
import logging

from . import _misc as _misc

"""
Access to network databases.
"""


def omnipath_universe(**kwargs):
    """
    Access generic networks from OmniPath.
    """

    return op.interactions.PostTranslational.get(**kwargs)


