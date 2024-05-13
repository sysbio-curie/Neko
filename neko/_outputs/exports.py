import pandas as pd


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
        net.convert_edgelist_into_genesymbol()
        self.nodes = net.nodes
        self.interactions = net.edges
        return

    def export_bnet(self, file_name="logic_model.bnet"):
        """
        Function to export the network in bnet format.
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
            # print source and target for each undefined interaction and relative references
            print("Bimodal interactions:")
            for index, row in undefined_interactions.iterrows():
                print(f"{row['source']} -> {row['target']}")
                print(f"Reference: {row['References']}")
        # Identify bimodal interactions
        bimodal_interactions = self.interactions.query("Effect == 'bimodal'")
        if not bimodal_interactions.empty:
            print(f"Warning: The network has {len(bimodal_interactions)} BIMODAL interaction(s).")
            # print source and target for each bimodal interaction and relative references
            print("Bimodal interactions:")
            for index, row in bimodal_interactions.iterrows():
                print(f"{row['source']} -> {row['target']}")
                print(f"Reference: {row['References']}")

        # Pre-filter stimulations, inhibitions, and exclude undefined effects
        stimulations = self.interactions.query("Effect == 'stimulation' or Effect == 'bimodal'")
        inhibitions = self.interactions.query("Effect == 'inhibition' or Effect == 'bimodal'")
        complex_formation = self.interactions.query("Effect == 'form complex'")

        with open(file_name, "w") as f:
            f.write("# model in BoolNet format\n")
            f.write("targets, factors\n")

            for entry in self.nodes.values:
                node = entry[0]
                formula_on = stimulations[stimulations["target"] == node]["source"].to_list()
                formula_off = inhibitions[inhibitions["target"] == node]["source"].to_list()
                formula_complex = complex_formation[complex_formation["target"] == node]["source"].to_list()

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

    def export_sif(self, file_name="logic_model.sif"):
        """
        Function to export the network in SIF format
        """

        with open(file_name, 'w') as file:
            for index, row in self.interactions.iterrows():
                # Use the Effect column directly assuming it contains "activate" or "inhibit"
                interaction_type = row['Effect']

                # Reference for the interaction
                interaction_reference = row['References']  # Adjust column name if necessary

                # Write a comment line with the interaction reference
                file.write(f"# Reference PMID: {interaction_reference}\n")

                # Write the formatted interaction to the .sif file
                file.write(f"{row['source']}\t{interaction_type}\t{row['target']}\n")

        return
