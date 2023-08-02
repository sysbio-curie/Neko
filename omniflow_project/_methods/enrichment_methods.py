import pandas as pd
from typing_extensions import Literal
from .._inputs.resources import Resources


class Connections:
    """
    Class that stores many utility functions to enrich an object Network.
    Each utility functions should take as input the nodes dataframe, which is used as base for each algorithm, and a
    database from the inputs modules, which will be used to extend the initial network.
    """

    def __init__(self):
        self.resources: Resources = Resources()

    def find_neighbours(self,
                        node=str,
                        mode: Literal['OUT', 'IN', 'ALL'] = 'ALL') -> list[str]:
        """
        helper function that finds the neighbours of the target node
        """

        db = self.resources.all_omnipath_interactions()

        targets = [i for i in db.loc[db["source"] == node][
            "target"]]  # storing all the target for the selected node
        sources = [j for j in db.loc[db["target"] == node][
            "source"]]  # storing all the target for the selected node
        if mode == "IN":
            return sources
        elif mode == "OUT":
            return targets
        else:
            return sources + targets

    def find_paths(self,
                   start=str,
                   end=str,
                   maxlen: int = 2,
                   minlen: int = 1,
                   loops: bool = False,
                   mode: Literal['OUT', 'IN', 'ALL'] = 'OUT',
                   ) -> list[tuple]:
        """
        Find paths or motifs in a network.
        Adapted from pypath function 'find_paths' in Network core class.
        In future this class will be extended to take into account many other parameters to select those paths that
        match certain criteria (type of interaction, effect, direction, data model, etc...)
        For now (version 0.1.0) it will act as skeleton to retrieve interactions from AllInteractions in Omnipath.
        """
        def find_all_paths_aux(start,
                               end,
                               path,
                               maxlen=None):

            path = path + [start]

            if (
                len(path) >= minlen + 1 and
                (
                    start == end or
                    (
                        end is None and
                        not loops and
                        len(path) == maxlen + 1
                    ) or
                    (
                        loops and
                        path[0] == path[-1]
                    )
                )
            ):
                return [path]

            paths = []

            if len(path) <= maxlen:

                next_steps = set(
                    self.find_neighbours(
                        start,
                        mode
                    )
                )

                next_steps = next_steps if loops else next_steps - set(path)

                for node in next_steps:
                    paths.extend(
                        find_all_paths_aux(
                            node,
                            end,
                            path, maxlen
                        )
                    )

            return paths

        minlen = max(1, minlen)

        all_paths = []

        all_paths.extend(find_all_paths_aux(start, end, [], maxlen))

        return all_paths
