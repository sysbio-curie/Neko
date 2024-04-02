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

    def find_neighbours(self, node: str, mode: Literal['OUT', 'IN']) -> list[str]:
        """
        Find the neighbours of a node in the network.
        Parameters:
        - node: Node for which to find neighbours.
        - mode: Direction of interactions to consider ('OUT' or 'IN').
        Returns:
        - neighbors: List of neighbouring nodes.
        """
        db = self.resources

        if mode == 'IN':
            neighbors = db.loc[db["target"] == node]["source"].tolist()
        else:
            neighbors = db.loc[db["source"] == node]["target"].tolist()

        return neighbors

    def find_paths(self,
                   start: (
                       str | pd.DataFrame | list[str]
                   ),
                   end: (
                       str | pd.DataFrame | list[str] | None
                   ) = None,
                   maxlen: int = 2,
                   minlen: int = 1,
                   loops: bool = False,
                   mode: Literal['OUT', 'IN'] = "OUT",
                   ) -> list[tuple]:
        """
        Find paths or motifs in a network.
        Adapted from pypath function 'find_paths' in Network core class.
        In future this class will be extended to take into account many other parameters to select those paths that
        match certain criteria (type of interaction, effect, direction, data model, etc...)
        For now (version 0.1.0) it will act as skeleton to retrieve interactions from AllInteractions in Omnipath.
        """

        def convert_to_string_list(start):
            """
            Convert the 'start' variable to a list of strings.
            """
            if isinstance(start, str):
                return [start]
            elif isinstance(start, pd.DataFrame):
                return start['name_of_node'].tolist()
            elif isinstance(start, list) and all(isinstance(item, str) for item in start):
                return start
            else:
                raise ValueError("Invalid type for 'start' variable")

        def find_all_paths_aux(start, end, path, maxlen):
            """
            Recursive function to find all paths between two nodes.
            Parameters:
            - start: Starting node.
            - end: Ending node.
            - path: Current path.
            - maxlen: Maximum length of paths to consider.
            Returns:
            - paths: List of all paths found.
            """
            path = path + [start]

            if len(path) >= minlen + 1 and (start == end or (end is None and not loops and len(path) == maxlen + 1) or (
                loops and path[0] == path[-1])):
                return [path]

            paths = []

            if len(path) <= maxlen:
                next_steps = self.find_neighbours(start, mode)

                if not loops:
                    next_steps = list(set(next_steps) - set(path))

                for node in next_steps:
                    paths.extend(find_all_paths_aux(node, end, path, maxlen))

            return paths

        start_nodes = convert_to_string_list(start)
        end_nodes = convert_to_string_list(end) if end else [None]

        minlen = max(1, minlen)
        all_paths = []

        for s in start_nodes:
            for e in end_nodes:
                all_paths.extend(find_all_paths_aux(s, e, [], maxlen))

        return all_paths
