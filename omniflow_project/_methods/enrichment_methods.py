from __future__ import annotations
from future.utils import iteritems
import pandas as pd
from typing_extensions import Literal


class Connections:
    """
    Class that stores many utility functions to enrich an object Network.
    Each utility functions should take as input the nodes dataframe, which is used as base for each algorithm, and a
    database from the inputs modules, which will be used to extend the initial network.
    """

    def __init__(self, database):
        self.resources = database
        return

    def find_neighbours(self,
                        node: str,
                        mode: Literal['OUT', 'IN', 'ALL'] = 'ALL') -> list[str]:
        """
        Optimized helper function that finds the neighbors of the target node.
        """

        db = self.resources

        target_nodes = set(db.loc[db["source"] == node]["target"])
        source_nodes = set(db.loc[db["target"] == node]["source"])

        if mode == "IN":
            return list(source_nodes)
        elif mode == "OUT":
            return list(target_nodes)
        else:
            return list(source_nodes | target_nodes)

    def find_paths(self,
                   start: (
                       str | pd.DataFrame | list[str]
                   ),
                   end: (
                       str | pd.DataFrame | None
                   ) = None,
                   maxlen: int = 2,
                   minlen: int = 1,
                   loops: bool = False,
                   mode: Literal['OUT', 'IN', 'ALL'] = 'ALL',
                   ) -> list[tuple]:
        """
        Find paths or motifs in a network.
        Adapted from pypath function 'find_paths' in Network core class.
        In future this class will be extended to take into account many other parameters to select those paths that
        match certain criteria (type of interaction, effect, direction, data model, etc...)
        For now (version 0.1.0) it will act as skeleton to retrieve interactions from AllInteractions in Omnipath.
        """

        def convert_to_string_list(start):
            if isinstance(start, str):
                return [start]
            elif isinstance(start, pd.DataFrame):
                return start['name_of_node'].tolist()
            elif isinstance(start, list) and all(isinstance(item, str) for item in start):
                return start
            else:
                raise ValueError("Invalid type for 'start' variable")

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

        start = convert_to_string_list(start)

        end = convert_to_string_list(end) if end else (None,)

        for s in start:

            for e in end:

                all_paths.extend(find_all_paths_aux(s, e, [], maxlen))

        return all_paths
