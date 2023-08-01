import omnipath as op
import pypath

class Resources():
    """
    This class stores the actual databases to mine interesting interactions. The user can select different
    database format from omnipath, pypath or loading one of its own.
    """

    def __init__(self):
        self.all_interactions = op.interactions.AllInteractions.get()

    def all_omnipath_interactions(self):
        return self.all_interactions
