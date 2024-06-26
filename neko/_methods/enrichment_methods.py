from typing import Union, List, Tuple
import pandas as pd
from collections import deque
import random


class Connections:
    """
    Class that stores many utility functions to enrich an object Network.
    Each utility functions should take as input the nodes dataframe, which is used as base for each algorithm, and a
    database from the inputs modules, which will be used to extend the initial network.
    """

    def __init__(self, database: pd.DataFrame):
        self.resources = database.copy()
        self.target_neighbours_map = self._preprocess_target_neighbours()
        self.source_neighbours_map = self._preprocess_source_neighbours()

    def _preprocess_target_neighbours(self) -> dict:
        """
        Preprocess the targets neighbours map for fast lookup.
        """
        df_sources = self.resources.set_index('source')
        return df_sources.groupby('source')['target'].apply(set).to_dict()

    def _preprocess_source_neighbours(self) -> dict:
        """
        Preprocess the sources neighbours map for fast lookup.
        """
        df_targets = self.resources.set_index('target')
        return df_targets.groupby('target')['source'].apply(set).to_dict()

    def find_target_neighbours(self, node: str) -> List[str]:
        """
        Optimized helper function that finds the neighbors of the target node.
        """
        return list(self.target_neighbours_map.get(node, []))

    def find_source_neighbours(self, node: str) -> List[str]:
        """
        Optimized helper function that finds the neighbors of the target node.
        """
        return list(self.source_neighbours_map.get(node, []))

    def find_all_neighbours(self, node: str) -> List[str]:
        """
        Optimized helper function that finds all neighbors (both source and target) of the target node.
        """
        target_neighs = self.find_target_neighbours(node)
        source_neighs = self.find_source_neighbours(node)
        return list(set(target_neighs + source_neighs))

    def find_paths(self,
                   start: Union[str, pd.DataFrame, List[str]],
                   end: Union[str, pd.DataFrame, List[str], None] = None,
                   maxlen: int = 2,
                   minlen: int = 1,
                   loops: bool = False) -> List[Tuple[str, ...]]:
        """
        Find paths or motifs in a network.
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

        def find_all_paths_aux(start, end, path, maxlen):
            path = path + [start]

            if len(path) >= minlen + 1 and (start == end or (end is None and not loops and len(path) == maxlen + 1) or (
                loops and path[0] == path[-1])):
                return [path]

            paths = []

            if len(path) <= maxlen:
                next_steps = self.find_target_neighbours(start)

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

    def bfs(self, start: str, end: str) -> List[Tuple[str, ...]]:
        """
        Find a path between two nodes using Breadth First Search. Useful to quickly check if two nodes are connected.

        Args:
            start: The start node.
            end: The end node.

        Returns:
            List of nodes representing the path from start to end. If no path exists, returns an empty list.
        """
        visited = set()
        queue = deque([(start, [start])])  # Store the path along with the node

        while queue:
            node, path = queue.popleft()
            if node == end:
                return [(path[i], path[i+1]) for i in range(len(path)-1)]

            if node not in visited:
                visited.add(node)
                neighbours = self.find_target_neighbours(node)
                random.shuffle(neighbours)  # Shuffle the neighbours to avoid bias
                queue.extend((neighbour, path + [neighbour]) for neighbour in neighbours if neighbour not in visited)

        return []

    def find_shortest_path(self,
                           start: Union[str, pd.DataFrame, List[str]],
                           end: Union[str, pd.DataFrame, List[str], None] = None,
                           loops: bool = False,
                           max_len: int = 6) -> List[Tuple[str, ...]]:
        """
        Find (one of) the shortest paths between two nodes in the network. It stops exploring path if it exceeds max_len.
        Slower than bfs function but ensures the shortest path.

        Args:
            start:
            end:
            loops:
            max_len: Maximum length of paths to explore.

        Returns:

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

        def find_shortest_path_aux(start, end, path, max_len):
            path = path + [start]

            if start == end or (end is None and not loops and len(path) == max_len + 1) or (
                loops and path[0] == path[-1]):
                return path

            shortest_path = []
            if len(path) <= max_len:
                next_steps = self.find_target_neighbours(start)

                if not loops:
                    next_steps = list(set(next_steps) - set(path))

                for node in next_steps:
                    newpath = find_shortest_path_aux(node, end, path, max_len)

                    if newpath:
                        if not shortest_path or len(newpath) < len(shortest_path): # here we discard the longest paths
                            shortest_path = newpath

            return shortest_path

        start_nodes = convert_to_string_list(start)
        end_nodes = convert_to_string_list(end) if end else [None]

        shortest_paths = []

        for s in start_nodes:
            for e in end_nodes:
                shortest_paths.extend(find_shortest_path_aux(s, e, [], max_len))

        return shortest_paths

    def find_all_shortest_paths(self,
                                start: Union[str, pd.DataFrame, List[str]],
                                end: Union[str, pd.DataFrame, List[str], None] = None,
                                loops: bool = False,
                                max_len: int = 6) -> List[Tuple[str, ...]]:
        """
        Find all shortest paths between two nodes in the network. It stops exploring path if it exceeds max_len.
        Slower than find_shortest_path function as it needs to keep track of all shortest paths.

        Args:
            start:
            end:
            loops:
            max_len: Maximum length of paths to explore.

        Returns:
            A list of all shortest paths.
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

        def find_all_shortest_paths_aux(start, end, path):
            path = path + [start]

            if start == end or (end is None and not loops and len(path) == max_len + 1) or (
                loops and path[0] == path[-1]):
                return [path]

            shortest_paths = []

            if len(path) <= max_len:
                next_steps = self.find_target_neighbours(start)

                if not loops:
                    next_steps = list(set(next_steps) - set(path))

                for node in next_steps:
                    newpaths = find_all_shortest_paths_aux(node, end, path)

                    for newpath in newpaths:
                        if not shortest_paths or len(newpath) == len(shortest_paths[0]):
                            shortest_paths.append(newpath)
                        elif len(newpath) < len(shortest_paths[0]):
                            shortest_paths = [newpath]

            return shortest_paths

        start_nodes = convert_to_string_list(start)
        end_nodes = convert_to_string_list(end) if end else [None]
        all_shortest_paths = []

        for s in start_nodes:
            for e in end_nodes:
                all_shortest_paths.extend(find_all_shortest_paths_aux(s, e, []))

        return all_shortest_paths

    def find_upstream_cascades(self,
                               target_genes: List[str],
                               max_depth: int = 1,
                               selected_rank: int = 1) -> List[Tuple[str, str]]:
        """
        Find cascades of interactions in the network.
        Parameters:
        - target_genes: List of target genes to start the cascade.
        - max_depth: Maximum depth of the cascade.
        - selected_rank: Number of top regulators to select for each iteration.
        Returns:
        - interactions: List of interactions in the cascade.
        """

        def collect_for_depth(current_targets, current_depth):
            if current_depth > max_depth:
                return []

            mcs_regulators = find_minimal_covering_regulators(self.resources, current_targets, selected_rank)

            interactions = [(reg, target) for reg in mcs_regulators for target in self.target_neighbours_map.get(reg, []) if target in current_targets]  # this is not working

            if current_depth < max_depth:
                next_targets = list(mcs_regulators)
                interactions += collect_for_depth(next_targets, current_depth + 1)

            return interactions

        return collect_for_depth(target_genes, 1) # it returns nothing


def remove_keys_and_related_values(dictionary: dict, keys_to_remove: list[str]) -> dict:
    """
    Remove keys and related values from a dictionary.
    """
    values_to_remove = set()
    for key in keys_to_remove:
        values_to_remove.update(dictionary.get(key, set))
        # Identify all keys to be removed based on the collected values
    all_keys_to_remove = set(keys_to_remove)
    for key, values in dictionary.items():
        if not values.isdisjoint(values_to_remove):  # If there's any overlap in values
            all_keys_to_remove.add(key)

    # Remove all identified keys from the dictionary
    for key in all_keys_to_remove:
        dictionary.pop(key, None)

    return dictionary


def check_dataframe_coverage(df: pd.DataFrame, target_genes: List[str]) -> List[str]:
    """
    Checks the DataFrame for the presence of the specified target genes and returns
    an updated list with those genes that are covered.
    """
    present_targets = set(df['target'].unique())
    targets = set(target_genes)
    covered_targets = list(targets & present_targets)
    return covered_targets


def find_minimal_covering_regulators(df: pd.DataFrame, target_genes: list[str], rank: int = 1):
    """
    Find the minimal set of regulators that cover all target genes in the DataFrame.

    Parameters:
    - df: pandas DataFrame with at least 'source' and 'target' columns.
    - target_genes: List of target genes to cover.
    - rank: Number of top regulators to select for each iteration.

    Returns:
    - regulators: Set of regulators that cover all target genes.
    """
    # Update target_genes list based on DataFrame coverage
    covered_targets = check_dataframe_coverage(df, target_genes)

    if len(covered_targets) < len(target_genes):
        missing_targets = set(target_genes) - set(covered_targets)
        print(f"Warning: Some target genes are not present in the DataFrame and will be ignored: {missing_targets}")

    # Proceed to filter the DataFrame for covered target genes
    filtered_df = df[df['target'].isin(covered_targets)]
    regulators = set()
    targets_covered = set()

    # Find unique regulators for each gene
    regulator_to_targets = filtered_df.groupby('source')['target'].apply(set).to_dict()

    while targets_covered != set(covered_targets):

        previous_covered_count = len(targets_covered)  # Keep track of previously covered targets
        # Calculate the lengths of values for each regulator
        lengths = sorted({len(values) for values in regulator_to_targets.values()}, reverse=True)

        # Select the top_n lengths
        selected_lengths = lengths[:rank]

        # Find all keys that have values lengths in the selected_lengths
        keys_with_max_values = [key for key, values in regulator_to_targets.items() if
                                len(values) in selected_lengths]

        regulators.update(keys_with_max_values)

        for reg in keys_with_max_values:
            targets_covered.update(regulator_to_targets[reg])

        # Check if no new targets were covered in this iteration
        if len(targets_covered) == previous_covered_count:
            print("Warning: Unable to cover all targets with available regulators. Exiting...")
            break
        regulator_to_targets = remove_keys_and_related_values(regulator_to_targets, keys_with_max_values)

    return regulators
