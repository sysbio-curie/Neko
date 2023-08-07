import omnipath as op
import pypath
from pypath.utils import mapping

class Resources():
    """
    This class stores the actual databases to mine interesting interactions. The user can select different
    database format from omnipath, pypath or loading one of its own.
    """

    def __init__(self):
        self.all_omnipath_interactions = None
        self.omnipath_interactions = None

    def load_all_omnipath_interactions(self):
        """
        loads into the Resources object the omnipath dataframe "all_interactions"
        """
        self.all_omnipath_interactions = op.interactions.AllInteractions.get()
        return

    def translate_dataframe_from_uniprot_to_genesymbol(self):
        """If the loaded dataframe has the source and target columns in uniprot, it will translate it to genesymbol"""

        return

    def load_omnipath_interactions(self):
        """
        loads into the Resources object the omnipath dataframe "Omnipath"
        """
        self.omnipath_interactions = op.interactions.OmniPath.get()
        return
