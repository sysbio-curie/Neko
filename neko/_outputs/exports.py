import os
import re
import itertools

import pandas as pd

from neko.core.tools import consolidate_edges


def _sanitize_bnet_identifier(value: str) -> str:
    """Convert a network label into a BoolNet-compatible identifier."""

    return re.sub(r"[\/\-\s\#:]", "_", str(value))


class Exports:
    """
    This class implement many methods used to export the Network object in different format.
    In particular the exports format will be methods-oriented (MaBoSS, Ginsim, cobrexa and so on...).
    To start with, the user can export the Network in SIF and Bnet format.
    In the future many more versatile methods will be implemented (SBML) and annotations will be included for each
    interaction, including the DOI of the relative reference and annotations from each database
    """
    def __init__(self, network):
        net = network.copy()
        df_edges = net.convert_edgelist_into_genesymbol()
        self.nodes = net.nodes
        self.interactions = consolidate_edges(df_edges)
        return

    def export_bnet(self, file_name="logic_model.bnet", n=None):
        """
        Function to export the network in bnet format, creating multiple files for bimodal interactions.
        """
        # Checks for nodes and interactions data
        if not isinstance(self.nodes, pd.DataFrame) or self.nodes.empty:
            print("Error: Nodes data is missing or empty.")
            return
        if not isinstance(self.interactions, pd.DataFrame) or self.interactions.empty:
            print("Error: Interactions data is missing or empty.")
            return

        required_node_columns = {'Genesymbol'}
        required_edge_columns = {'source', 'target', 'Effect'}
        missing_node_columns = required_node_columns.difference(
            self.nodes.columns,
        )
        missing_edge_columns = required_edge_columns.difference(
            self.interactions.columns,
        )
        if missing_node_columns or missing_edge_columns:
            missing = sorted(missing_node_columns | missing_edge_columns)
            raise ValueError(
                'BNet export is missing required column(s): '
                f'{", ".join(missing)}.',
            )

        interactions = consolidate_edges(self.interactions)

        invalid_nodes = self.nodes['Genesymbol'].isna() | self.nodes[
            'Genesymbol'
        ].map(lambda value: isinstance(value, str) and not value.strip())
        invalid_sources = interactions['source'].isna() | interactions[
            'source'
        ].map(lambda value: isinstance(value, str) and not value.strip())
        invalid_targets = interactions['target'].isna() | interactions[
            'target'
        ].map(lambda value: isinstance(value, str) and not value.strip())

        if invalid_nodes.any() or invalid_sources.any() or invalid_targets.any():
            raise ValueError(
                'BNet export requires non-empty node labels and interaction '
                'endpoints; found a null or empty identifier.',
            )

        node_labels = list(dict.fromkeys(self.nodes['Genesymbol'].tolist()))
        node_set = set(node_labels)
        edge_nodes = set(interactions['source']).union(interactions['target'])
        missing_nodes = edge_nodes.difference(node_set)
        if missing_nodes:
            detail = ', '.join(sorted(map(str, missing_nodes)))
            raise ValueError(
                'BNet interaction endpoints are absent from the node table: '
                f'{detail}.',
            )

        # Identify undefined interactions
        undefined_interactions = interactions.query("Effect == 'undefined'")
        if not undefined_interactions.empty:
            print(f"Warning: The network has {len(undefined_interactions)} UNDEFINED interaction(s).")
            print("Undefined interactions:")
            for _, row in undefined_interactions.iterrows():
                print(f"{row['source']} -> {row['target']}")
                print(f"Reference: {row['References']}")

        # Identify bimodal interactions
        bimodal_interactions = interactions.query("Effect == 'bimodal'")
        if not bimodal_interactions.empty:
            print(f"Warning: The network has {len(bimodal_interactions)} BIMODAL interaction(s).")
            print("Bimodal interactions:")
            for _, row in bimodal_interactions.iterrows():
                print(f"{row['source']} -> {row['target']}")
                print(f"Reference: {row['References']}")

        if n is not None and (not isinstance(n, int) or isinstance(n, bool) or n < 0):
            raise ValueError('n must be a non-negative integer or None.')

        # Keep this iterator lazy: materializing all 2^k variants defeats the
        # purpose of `n` and can exhaust memory before the first file is made.
        permutations = itertools.product(
            ['stimulation', 'inhibition'],
            repeat=len(bimodal_interactions),
        )

        if n is not None:
            permutations = itertools.islice(permutations, n)

        bimodal_indices = bimodal_interactions.index.tolist()

        sanitized_nodes = {
            node: _sanitize_bnet_identifier(node)
            for node in node_labels
        }
        collisions = {}

        for original, sanitized in sanitized_nodes.items():
            collisions.setdefault(sanitized, []).append(original)

        collisions = {
            sanitized: originals
            for sanitized, originals in collisions.items()
            if len(originals) > 1
        }

        if collisions:
            detail = '; '.join(
                f'{sanitized}: {", ".join(map(str, originals))}'
                for sanitized, originals in collisions.items()
            )
            raise ValueError(
                'Node labels collide after BNet identifier sanitization '
                f'({detail}).',
            )

        # Create a directory for the BNet files if a directory is provided
        directory = os.path.dirname(file_name)
        if directory:
            os.makedirs(directory, exist_ok=True)

        # Iterate through permutations and create a BNet file for each
        generated = 0

        for i, perm in enumerate(permutations):
            # Create a copy of the interactions DataFrame
            interactions_copy = interactions.copy()

            # Update bimodal interactions based on the current permutation
            for interaction_index, effect in zip(bimodal_indices, perm):
                interactions_copy.loc[interaction_index, 'Effect'] = effect

            # Pre-filter stimulations, inhibitions, and exclude undefined effects
            stimulations = interactions_copy.query("Effect == 'stimulation'")
            inhibitions = interactions_copy.query("Effect == 'inhibition'")
            complex_formation = interactions_copy.query("Effect == 'form complex'")

            # Generate the file name for this permutation
            perm_file_name = f"{os.path.splitext(file_name)[0]}_{i + 1}.bnet"

            with open(perm_file_name, "w") as f:
                f.write("# model in BoolNet format\n")
                f.write("targets, factors\n")

                for node in node_labels:
                    sanitized_node = sanitized_nodes[node]
                    formula_on = [
                        _sanitize_bnet_identifier(src)
                        for src in stimulations.loc[
                            stimulations['target'] == node,
                            'source',
                        ].tolist()
                    ]
                    formula_off = [
                        _sanitize_bnet_identifier(src)
                        for src in inhibitions.loc[
                            inhibitions['target'] == node,
                            'source',
                        ].tolist()
                    ]
                    formula_complex = [
                        _sanitize_bnet_identifier(src)
                        for src in complex_formation.loc[
                            complex_formation['target'] == node,
                            'source',
                        ].tolist()
                    ]

                    # Constructing the formula
                    formula_parts = []
                    if formula_complex:
                        formula_parts.append(f"({' & '.join(formula_complex)})")
                    if formula_on:
                        formula_parts.append(f"({' | '.join(formula_on)})")
                    if formula_off:
                        formula_parts.append("!({})".format(" | ".join(formula_off)))

                    # Writing the node and its formula to the file
                    formula = ' & '.join(formula_parts) if formula_parts else sanitized_node
                    f.write(f"{sanitized_node}, {formula}\n")

            print(f"Created BNet file: {perm_file_name}")
            generated += 1

        print(f"Generated {generated} BNet files.")

    def export_sif(self, file_name="logic_model.sif"):
        """
        Function to export the network in SIF format
        """

        directory = os.path.dirname(file_name)

        if directory:
            os.makedirs(directory, exist_ok=True)

        with open(file_name, 'w') as file:
            for _, row in self.interactions.iterrows():
                # Use the Effect column directly assuming it contains "activate" or "inhibit"
                interaction_type = row['Effect']

                if interaction_type == "form complex":
                    interaction_type = "form_complex"

                # Reference for the interaction
                interaction_reference = row['References']  # Adjust column name if necessary

                # Write a comment line with the interaction reference
                file.write(f"# Reference PMID: {interaction_reference}\n")

                # Write the formatted interaction to the .sif file
                file.write(f"{row['source']}\t{interaction_type}\t{row['target']}\n")

        return
