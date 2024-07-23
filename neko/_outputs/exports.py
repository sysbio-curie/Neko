import pandas as pd
import os
import itertools
import re

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
        self.interactions = df_edges
        return

    def export_bnet(self, file_name="logic_model.bnet"):
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

        # Identify undefined interactions
        undefined_interactions = self.interactions.query("Effect == 'undefined'")
        if not undefined_interactions.empty:
            print(f"Warning: The network has {len(undefined_interactions)} UNDEFINED interaction(s).")
            print("Undefined interactions:")
            for index, row in undefined_interactions.iterrows():
                print(f"{row['source']} -> {row['target']}")
                print(f"Reference: {row['References']}")

        # Identify bimodal interactions
        bimodal_interactions = self.interactions.query("Effect == 'bimodal'")
        if not bimodal_interactions.empty:
            print(f"Warning: The network has {len(bimodal_interactions)} BIMODAL interaction(s).")
            print("Bimodal interactions:")
            for index, row in bimodal_interactions.iterrows():
                print(f"{row['source']} -> {row['target']}")
                print(f"Reference: {row['References']}")

        # Generate permutations for bimodal interactions
        bimodal_sources = bimodal_interactions['source'].tolist()
        bimodal_targets = bimodal_interactions['target'].tolist()
        permutations = list(itertools.product(['stimulation', 'inhibition'], repeat=len(bimodal_interactions)))

        # Create a directory for the BNet files
        os.makedirs(os.path.dirname(file_name), exist_ok=True)

        # Iterate through permutations and create a BNet file for each
        for i, perm in enumerate(permutations):
            # Create a copy of the interactions DataFrame
            interactions_copy = self.interactions.copy()

            # Update bimodal interactions based on the current permutation
            for j, (source, target) in enumerate(zip(bimodal_sources, bimodal_targets)):
                interactions_copy.loc[(interactions_copy['source'] == source) &
                                      (interactions_copy['target'] == target), 'Effect'] = perm[j]

            # Pre-filter stimulations, inhibitions, and exclude undefined effects
            stimulations = interactions_copy.query("Effect == 'stimulation'")
            inhibitions = interactions_copy.query("Effect == 'inhibition'")
            complex_formation = interactions_copy.query("Effect == 'form complex'")

            # Generate the file name for this permutation
            perm_file_name = f"{os.path.splitext(file_name)[0]}_{i + 1}.bnet"

            with open(perm_file_name, "w") as f:
                f.write("# model in BoolNet format\n")
                f.write("targets, factors\n")

                for entry in self.nodes.values:
                    node = entry[0]

                    # Replace special characters in node names
                    node = re.sub(r"[\/\-\s\#]", "_", node)

                    formula_on = [re.sub(r"[\/\-\s\#]", "_", src) for src in
                                  stimulations[stimulations["target"] == node]["source"].to_list()]
                    formula_off = [re.sub(r"[\/\-\s\#]", "_", src) for src in
                                   inhibitions[inhibitions["target"] == node]["source"].to_list()]
                    formula_complex = [re.sub(r"[\/\-\s\#]", "_", src) for src in
                                       complex_formation[complex_formation["target"] == node]["source"].to_list()]

                    # Constructing the formula
                    formula_parts = []
                    if formula_complex:
                        formula_parts.append(f"({' & '.join(formula_complex)})")
                    if formula_on:
                        formula_parts.append(f"({' | '.join(formula_on)})")
                    if formula_off:
                        formula_parts.append("!({})".format(" | ".join(formula_off)))

                    # Writing the node and its formula to the file
                    f.write(f"{node}, {' & '.join(formula_parts) if formula_parts else node}\n")

            print(f"Created BNet file: {perm_file_name}")

        print(f"Generated {len(permutations)} BNet files.")

    def export_sif(self, file_name="logic_model.sif"):
        """
        Function to export the network in SIF format
        """

        with open(file_name, 'w') as file:
            for index, row in self.interactions.iterrows():
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
