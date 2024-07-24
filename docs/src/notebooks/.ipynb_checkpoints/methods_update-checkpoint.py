from typing import Union, List, Tuple
import pandas as pd


class Connections:
    """
    Class that stores many utility functions to enrich an object Network.
    Each utility functions should take as input the nodes dataframe, which is used as base for each algorithm, and a
    database from the inputs modules, which will be used to extend the initial network.
    """

    def __init__(self, database: pd.DataFrame):
        self.resources = database.copy()
        self.resources.set_index('source', inplace=True)
        self.neighbours_map = self._preprocess_neighbours()

    def _preprocess_neighbours(self) -> dict:
        """
        Preprocess the neighbours map for fast lookup.
        """
        return self.resources.groupby('source')['target'].apply(set).to_dict()

    def find_neighbours(self, node: str) -> List[str]:
        """
        Optimized helper function that finds the neighbors of the target node.
        """
        return list(self.neighbours_map.get(node, []))

    def find_all_neighbours(self, node: str) -> List[str]:
        """
        Optimized helper function that finds all neighbors (both source and target) of the target node.
        """
        target_neighs = self.resources.loc[self.resources.index == node, 'target'].tolist()
        source_neighs = self.resources[self.resources['target'] == node].index.tolist()
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

            if len(path) >= minlen + 1 and (start == end or (end is None and not loops and len(path) == maxlen + 1) or (loops and path[0] == path[-1])):
                return [path]

            paths = []

            if len(path) <= maxlen:
                next_steps = self.find_neighbours(start)

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

            interactions = [(reg, target) for reg in mcs_regulators for target in self.resources[self.resources.index == reg]['target'] if target in current_targets]

            if current_depth < max_depth:
                next_targets = list(mcs_regulators)
                interactions += collect_for_depth(next_targets, current_depth + 1)

            return interactions

        return collect_for_depth(target_genes, 1)


def remove_keys_and_related_values(dictionary: dict, keys_to_remove: List[str]) -> dict:
    """
    Remove keys and related values from a dictionary.
    """
    values_to_remove = set()
    for key in keys_to_remove:
        values_to_remove.update(dictionary.get(key, set()))

    all_keys_to_remove = set(keys_to_remove)
    for key, values in dictionary.items():
        if not values.isdisjoint(values_to_remove):
            all_keys_to_remove.add(key)

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


def find_minimal_covering_regulators(df: pd.DataFrame, target_genes: List[str], rank: int = 1) -> set:
    """
    Find the minimal set of regulators that cover all target genes in the DataFrame.
    """
    covered_targets = check_dataframe_coverage(df, target_genes)

    if len(covered_targets) < len(target_genes):
        missing_targets = set(target_genes) - set(covered_targets)
        print(f"Warning: Some target genes are not present in the DataFrame and will be ignored: {missing_targets}")

    filtered_df = df[df['target'].isin(covered_targets)]
    regulators = set()
    targets_covered = set()

    regulator_to_targets = filtered_df.groupby(filtered_df.index)['target'].apply(set).to_dict()

    while targets_covered != set(covered_targets):
        previous_covered_count = len(targets_covered)
        lengths = sorted({len(values) for values in regulator_to_targets.values()}, reverse=True)
        selected_lengths = lengths[:rank]
        keys_with_max_values = [key for key, values in regulator_to_targets.items() if len(values) in selected_lengths]
        regulators.update(keys_with_max_values)

        for reg in keys_with_max_values:
            targets_covered.update(regulator_to_targets[reg])

        if len(targets_covered) == previous_covered_count:
            print("Warning: Unable to cover all targets with available regulators. Exiting...")
            break

        regulator_to_targets = remove_keys_and_related_values(regulator_to_targets, keys_with_max_values)

    return regulators
