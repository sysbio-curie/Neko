from __future__ import annotations
import pandas as pd
from typing_extensions import Literal


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


def check_dataframe_coverage(df: pd.DataFrame, target_genes: list[str]) -> list[str]:
    """
    Checks the DataFrame for the presence of the specified target genes and returns
    an updated list with those genes that are covered.

    Parameters:
    - df: pandas DataFrame with at least 'target' column.
    - target_genes: Iterable of target genes to check for in the DataFrame.

    Returns:
    - covered_targets: List of target genes that are present in the DataFrame.
    """
    # Determine which target genes are present in the DataFrame
    present_targets = set(df['target'].unique())
    targets = set(target_genes)

    # Filter out target genes not present in the DataFrame
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
                        node: str) -> list[str]:
        """
        Optimized helper function that finds the neighbors of the target node.
        """

        db = self.resources

        neigh = set(db.loc[db["source"] == node]["target"])

        return list(neigh)

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

        def find_all_paths_aux(start, end, path, maxlen):
            path = path + [start]

            if len(path) >= minlen + 1 and (start == end or (end is None and not loops and len(path) == maxlen + 1) or (
                loops and path[0] == path[-1])):
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
                               target_genes: list[str],
                               max_depth: int = 1,
                               selected_rank: int = 1) -> list[tuple]:

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
            # Base case: Exit if maximum depth is exceeded
            if current_depth > max_depth:
                return []

            # Find the minimal covering set of regulators for current targets
            mcs_regulators = find_minimal_covering_regulators(self.resources, current_targets, selected_rank)

            # Collect interactions (source, target) for these regulators
            interactions = [(reg, target) for reg in mcs_regulators for target in self.resources[self.resources['source'] == reg]['target'] if
                            target in current_targets]

            # If additional depth is required, treat the regulators just found as the next targets
            if current_depth < max_depth:
                next_targets = mcs_regulators  # The regulators themselves become the targets for the next depth
                interactions += collect_for_depth(next_targets, current_depth + 1)

            return interactions

        # Kick off the recursive process starting with the original genes as targets
        return collect_for_depth(target_genes, 1)
