from __future__ import annotations
from pypath.utils import mapping
import pandas as pd
import networkx as nx

def is_connected(network) -> bool:
    """
    This function checks if a network is connected. It takes a Network object as input and returns True if the network
    is connected, otherwise it returns False.

    Args:
        - network: A Network object representing the network to be checked.

    Returns:
        - bool
    """
    # Create a graph from the edges
    g = nx.from_pandas_edgelist(network.edges, 'source', 'target')
    # Add isolated nodes to the graph
    all_nodes = set(network.nodes['Uniprot'])
    g.add_nodes_from(all_nodes)
    # Check if the graph is connected
    return nx.is_connected(g)

def check_sign(interaction: pd.DataFrame, consensus: bool = False) -> str:
    """
    This function checks the sign of an interaction in the Omnipath format (Pandas DataFrame or Series).
    The attribute "consensus" checks for the consistency of the sign of the interaction among the references.

    Args:
        - interaction: A pandas DataFrame or Series representing the interaction.
        - consensus: A boolean indicating whether to check for consensus among references.

    Returns:
        - A string indicating the sign of the interaction: "stimulation", "inhibition", "form complex", or "undefined".
    """
    # Handle both DataFrame and Series input
    if isinstance(interaction, pd.DataFrame):
        interaction = interaction.iloc[0]

    if consensus:
        if interaction.get("consensus_inhibition") and interaction.get("consensus_stimulation"):
            return "bimodal"
        if interaction.get("consensus_stimulation"):
            return "stimulation"
        elif interaction.get("consensus_inhibition"):
            return "inhibition"
        else:
            return "undefined"
    else:
        # Check if it is both stimulation and inhibition
        if interaction.get("is_stimulation", True) and interaction.get("is_inhibition", True):
            return "bimodal"
        if interaction.get("is_stimulation", False):
            return "stimulation"
        elif interaction.get("is_inhibition", False):
            return "inhibition"
        # Check for "form_complex" column existence
        elif interaction.get("form_complex", False):
            return "form complex"
        else:
            return "undefined"


def check_gene_list_format(gene_list: list[str]) -> bool:
    """
    This function checks the format of the gene list and returns True if the gene list is in Uniprot format,
    False if the gene list is in genesymbol format.

    Args:
        - gene_list: A list of gene identifiers. The gene identifiers can be either Uniprot identifiers or genesymbols.

    Returns:
        - A boolean indicating whether the gene list is in Uniprot format (True) or genesymbol format (False).
    """
    # Check if the gene list contains Uniprot identifiers
    if all(mapping.id_from_label0(gene) for gene in gene_list):
        return True
    # Check if the gene list contains genesymbols
    elif all(mapping.label(gene) for gene in gene_list):
        return False


def mapping_node_identifier(node: str) -> list[str]:
    """
    This function takes a node identifier and returns a list containing the possible identifiers for the node.
    The identifiers include a complex string, a genesymbol, and a uniprot identifier. The function uses the
    mapping.id_from_label0 and mapping.label functions from the pypath.utils.mapping module to translate the node
    identifier into these different formats.

    Args:
        - node: A string representing the node identifier. The node identifier can be a genesymbol, a uniprot identifier,
                or a complex string.

    Returns:
        - A list containing the complex string, genesymbol, and uniprot identifier for the node. If the node identifier
          cannot be translated into one of these formats, the corresponding value in the list is None.
    """
    complex_string = None
    genesymbol = None
    uniprot = None

    if mapping.id_from_label0(node):
        # Convert UniProt ID to gene symbol
        uniprot = mapping.id_from_label0(node)
        if uniprot.startswith("MI"):
            genesymbol = uniprot
        else:
            # Set the UniProt ID as the 'Uniprot' value in the new entry
            genesymbol = mapping.label(uniprot)
    elif mapping.id_from_label0(node).startswith("COMPLEX"):
        node = node[8:]
        node_list = node.split("_")

        # Translate each element in node_list using mapping.label
        translated_node_list = [mapping.label(mapping.id_from_label0(item)) for item in node_list]

        # Join the elements in node_list with "_"
        joined_node_string = "_".join(translated_node_list)

        # Add back the "COMPLEX:" prefix to the string
        complex_string = "COMPLEX:" + joined_node_string
    elif mapping.label(node):
        genesymbol = mapping.label(node)
        uniprot = mapping.id_from_label0(genesymbol)
    else:
        print("Error during translation, check syntax for ", node)

    return [complex_string, genesymbol, uniprot]


def translate_paths(paths) -> list[list[str]]:
    """
    This function translates a list of paths, where each path is a sequence of node identifiers.
    It uses the helper function `handle_complex_identifier` to translate each node identifier in the paths.

    Args:
        - paths: A list of paths, where each path is a sequence of node identifiers.
                 A node identifier can be a string or a list of strings.

    Returns:
        - A list of translated paths, where each path is a sequence of translated node identifiers.
    """
    translated_list = []

    def handle_complex_identifier(item):
        """
        This helper function translates a node identifier using the `mapping_node_identifier` function.
        It checks all possible identifiers (complex, genesymbol, uniprot) and returns the first non-None value.

        Args:
        - item: A node identifier.

        Returns:
        - The translated node identifier.
        """
        identifiers = mapping_node_identifier(item)
        return identifiers[0] or identifiers[1] or identifiers[2]

    # If input_list is a list of strings
    if isinstance(paths[0], str):
        translated_list = [handle_complex_identifier(item) for item in paths]
    # If input_list is a list of lists of strings
    elif isinstance(paths[0], list):
        for sublist in paths:
            translated_sublist = [handle_complex_identifier(item) for item in sublist]
            translated_list.append(translated_sublist)

    return translated_list


def join_unique(series) -> str:
    """
    This function takes a pandas Series, filters out None values, and returns a string of unique values joined by a comma.

    Args:
        - series: A pandas Series object.

    Returns: - A string of unique values in the series, joined by a comma. If a value in the series is None,
                it is not included in the output string.
    """
    # Filter out None values before converting to set and joining
    filtered_series = [str(item) for item in series if item is not None]
    unique_items = set(filtered_series)
    return ', '.join(unique_items)


def determine_most_frequent_effect(effects):
    effect_counts = effects.value_counts()

    # If there's only one type of effect, return it
    if len(effect_counts) == 1:
        return effect_counts.index[0]

    # Count the occurrences of each effect type
    stimulation_count = effect_counts.get('stimulation', 0)
    inhibition_count = effect_counts.get('inhibition', 0)
    form_complex_count = effect_counts.get('form_complex', 0)

    # Calculate the total of stimulation, inhibition, and form_complex
    total_known_effects = stimulation_count + inhibition_count + form_complex_count

    # If there are more unknown effects than known ones, return 'undefined'
    if len(effects) - total_known_effects > total_known_effects:
        return 'undefined'

    # If form_complex is the most frequent, return it
    if form_complex_count > stimulation_count and form_complex_count > inhibition_count:
        return 'form_complex'

    # Handle stimulation and inhibition
    if stimulation_count > inhibition_count:
        return 'stimulation'
    elif inhibition_count > stimulation_count:
        return 'inhibition'
    elif stimulation_count == inhibition_count and stimulation_count > 0:
        return 'bimodal'

    # If we've reached this point, it means there's no clear majority
    return 'undefined'
