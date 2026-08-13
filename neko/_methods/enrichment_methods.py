from typing import Union, List, Tuple,Optional
import pandas as pd
from collections import deque, defaultdict
import random
from neko.core.tools import check_sign


class Connections:
    """
    Class that stores many utility functions to enrich an object Network.
    Each utility functions should take as input the nodes dataframe, which is used as base for each algorithm, and a
    database from the inputs modules, which will be used to extend the initial network.
    """

    def __init__(self, database: pd.DataFrame):
        self.resources = database.copy()
        self.resource_nodes = set(self.resources['source'].dropna())
        self.resource_nodes.update(self.resources['target'].dropna())
        self._edge_positions = defaultdict(list)
        for position, (source, target) in enumerate(
                self.resources[['source', 'target']].itertuples(
                    index=False,
                    name=None,
                )):
            self._edge_positions[(source, target)].append(position)
        self.target_neighbours_map = self._preprocess_target_neighbours()
        self.source_neighbours_map = self._preprocess_source_neighbours()
        # Precompute edge sign cache for both consensus True/False. Resource
        # tables carry boolean evidence columns, while working network edges
        # carry one canonical ``Effect`` value.
        self.signed_edges = {}
        self.signed_edges_consensus = {}
        signed, signed_consensus = self._signed_masks()
        edge_keys = self.resources[['source', 'target']].itertuples(
            index=False,
            name=None,
        )
        for key, is_signed, is_consensus_signed in zip(
                edge_keys,
                signed,
                signed_consensus,
            ):
            # Preserve the historical last-row-wins behavior for parallel
            # resource rows. Resource adapters normally consolidate these.
            self.signed_edges[key] = bool(is_signed)
            self.signed_edges_consensus[key] = bool(is_consensus_signed)

    def _signed_masks(self):
        resource_columns = {
            'is_stimulation',
            'is_inhibition',
            'form_complex',
        }
        if resource_columns.issubset(self.resources.columns):
            stimulation = self.resources['is_stimulation'].fillna(False).astype(bool)
            inhibition = self.resources['is_inhibition'].fillna(False).astype(bool)
            complex_formation = self.resources['form_complex'].fillna(False).astype(bool)
            signed = stimulation | inhibition | complex_formation

            if {
                    'consensus_stimulation',
                    'consensus_inhibition',
                }.issubset(self.resources.columns):
                signed_consensus = (
                    self.resources['consensus_stimulation'].fillna(False).astype(bool)
                    | self.resources['consensus_inhibition'].fillna(False).astype(bool)
                )
            else:
                signed_consensus = pd.Series(
                    False,
                    index=self.resources.index,
                )
            return signed, signed_consensus

        effects = self.resources.get(
            'Effect',
            pd.Series('undefined', index=self.resources.index),
        )
        effects = (
            effects
            .fillna('undefined')
            .astype(str)
            .str.strip()
            .str.lower()
            .str.replace('_', ' ', regex=False)
        )
        signed = effects.isin({
            'stimulation',
            'inhibition',
            'bimodal',
            'form complex',
        })
        return signed, signed.copy()

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
        return sorted(self.target_neighbours_map.get(node, []), key=str)

    def find_source_neighbours(self, node: str) -> List[str]:
        """
        Optimized helper function that finds the neighbors of the target node.
        """
        return sorted(self.source_neighbours_map.get(node, []), key=str)

    def find_all_neighbours(self, node: str) -> List[str]:
        """
        Optimized helper function that finds all neighbors (both source and target) of the target node.
        """
        target_neighs = self.find_target_neighbours(node)
        source_neighs = self.find_source_neighbours(node)
        return sorted(set(target_neighs + source_neighs), key=str)

    def find_interactions(self, source: str, target: str) -> pd.DataFrame:
        """Return resource rows for one directed edge without scanning the table."""

        positions = self._edge_positions.get((source, target), ())
        if not positions:
            return self.resources.iloc[0:0]
        return self.resources.iloc[positions]

    def is_signed_edge(self, source, target, consensus=False):
        """
        Returns True if the edge from source to target is signed (not undefined), False otherwise.
        Uses precomputed cache for speed.
        """
        key = (source, target)
        if consensus:
            result = self.signed_edges_consensus.get(key, False)
        else:
            result = self.signed_edges.get(key, False)
        return result

    def bfs(self, start: str, end: str, maxlen: Optional[int], only_signed: bool = False, consensus: bool = False, force: bool = False) -> List[List[str]]:
        """
        Returns the shortest path between two nodes (as a list of nodes) using BFS,
        but stops searching if the path length exceeds `maxlen` edges (if provided).
        If only_signed is True, only considers signed edges (not undefined).
        If force is False and maxlen is None, uses a default upper bound of 10.
        """
        if start == end:
            return [[start]]  # trivial path

        visited = set()
        # Set a default upper bound if maxlen is None and not force
        effective_maxlen = maxlen
        if maxlen is None and not force:
            effective_maxlen = 10
        queue = deque([(start, [start], 0)])  # (node, path_so_far, depth)

        while queue:
            node, path, depth = queue.popleft()
            if node == end:
                return [path]
            if effective_maxlen is not None and depth >= effective_maxlen:
                continue
            if node not in visited:
                visited.add(node)
                for neighbor in self.find_target_neighbours(node):
                    if neighbor not in visited:
                        if not only_signed or self.is_signed_edge(node, neighbor, consensus):
                            queue.append((neighbor, path + [neighbor], depth + 1))
        return []

    def bfs_all_shortest_edges(
            self,
            start: str,
            end: str,
            maxlen: int,
            only_signed: bool = False,
            consensus: bool = False,
        ) -> List[Tuple[str, str]]:
        """Return the edge union of all minimum-length paths within a cutoff."""

        if start == end:
            return []

        distances = {start: 0}
        predecessors = defaultdict(set)
        queue = deque([start])
        shortest_distance = None

        while queue:
            node = queue.popleft()
            depth = distances[node]

            if depth >= maxlen:
                continue
            if shortest_distance is not None and depth >= shortest_distance:
                continue

            for neighbor in self.find_target_neighbours(node):
                if only_signed and not self.is_signed_edge(
                        node,
                        neighbor,
                        consensus,
                    ):
                    continue

                next_depth = depth + 1
                known_depth = distances.get(neighbor)
                if known_depth is None:
                    distances[neighbor] = next_depth
                    predecessors[neighbor].add(node)
                    queue.append(neighbor)
                elif known_depth == next_depth:
                    predecessors[neighbor].add(node)

                if neighbor == end:
                    shortest_distance = next_depth

        if end not in distances:
            return []

        shortest_edges = set()
        pending = [end]
        expanded = set()
        while pending:
            node = pending.pop()
            if node in expanded:
                continue
            expanded.add(node)
            for predecessor in predecessors[node]:
                shortest_edges.add((predecessor, node))
                pending.append(predecessor)

        return sorted(shortest_edges, key=lambda edge: (str(edge[0]), str(edge[1])))


    def find_paths(self,
                   start: Union[str, pd.DataFrame, List[str]],
                   end: Union[str, pd.DataFrame, List[str], None] = None,
                   maxlen: int = 2,
                   minlen: int = 1,
                   loops: bool = False,
                   only_signed: bool = False,
                   consensus: bool = False) -> List[List[str]]:
        """
        Find all paths or motifs in a network, with optional sign/consensus filtering.
        Uses an iterative DFS with an explicit stack for better performance and memory efficiency.
        Args:
            start: Node(s) to start from (str, list of str, or DataFrame with 'name_of_node').
            end: Node(s) to end at (str, list of str, DataFrame, or None for motif search).
            maxlen: Maximum path length (number of edges).
            minlen: Minimum path length (number of edges).
            loops: Allow cycles/loops if True.
            only_signed: If True, only consider signed edges (not undefined).
            consensus: If True, use consensus sign filtering.
        Returns:
            List of paths (each path is a list of node names).
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

        def path_generator(start_nodes, end_nodes, maxlen, minlen, loops, only_signed, consensus):
            for s in start_nodes:
                for e in end_nodes:
                    stack = [(s, [s])]
                    while stack:
                        current, path = stack.pop()
                        # Prune if path too long
                        if len(path) > maxlen + 1:
                            continue
                        # Check for valid path
                        if len(path) >= minlen + 1 and (
                            (e is not None and current == e) or
                            (e is None and not loops and len(path) == maxlen + 1) or
                            (loops and path[0] == path[-1] and len(path) > 1)
                        ):
                            yield path
                        # Continue DFS
                        if len(path) <= maxlen:
                            next_steps = self.find_target_neighbours(current)
                            if only_signed:
                                next_steps = [n for n in next_steps if self.is_signed_edge(current, n, consensus)]
                            if not loops:
                                next_steps = list(set(next_steps) - set(path))
                            for neighbor in next_steps:
                                stack.append((neighbor, path + [neighbor]))

        start_nodes = convert_to_string_list(start)
        end_nodes = convert_to_string_list(end) if end else [None]
        minlen = max(1, minlen)
        # Collect all paths in a list for backward compatibility
        return list(path_generator(start_nodes, end_nodes, maxlen, minlen, loops, only_signed, consensus))

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
