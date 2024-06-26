from __future__ import annotations
from typing import List, Optional, Tuple
from pypath.utils import mapping
import omnipath as op
from itertools import combinations
import pandas as pd
from .._inputs.resources import Resources
from .._methods.enrichment_methods import Connections
from typing_extensions import Literal
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
import copy
from .._annotations.gene_ontology import Ontology

import pandas as pd

import networkx as nx


def is_connected(network) -> bool:
    """
    Return True if all the nodes in the nodes list in the Network object are connected, otherwise it returns False.
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

    Parameters:
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


def check_gene_list_format(gene_list: list[str]) -> list[str]:
    """
    This function checks the format of the gene list and returns True if the gene list is in Uniprot format, False if the gene list is in genesymbol format.

    Parameters:
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

    Parameters:
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


def translate_paths(paths):
    """
    This function translates a list of paths, where each path is a sequence of node identifiers.
    It uses the helper function `handle_complex_identifier` to translate each node identifier in the paths.

    Parameters:
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

        Parameters:
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


def join_unique(series):
    """
    This function takes a pandas Series, filters out None values, and returns a string of unique values joined by a comma.

    Parameters:
    - series: A pandas Series object.

    Returns:
    - A string of unique values in the series, joined by a comma. If a value in the series is None, it is not included in the output string.
    """
    # Filter out None values before converting to set and joining
    filtered_series = [str(item) for item in series if item is not None]
    unique_items = set(filtered_series)
    return ', '.join(unique_items)


class Network:
    """
    The Network class is the main class of the Omniflow package. It is designed to store nodes and edges of a network and offers various methods for enrichment analysis. The class takes a list of nodes as the main argument and a series of optional filters and other options such as inputs/outputs, specific tissues, type of interactions, prune inputs/outputs, etc.

    Attributes:
        nodes (pd.DataFrame): A DataFrame storing the nodes in the network.
        edges (pd.DataFrame): A DataFrame storing the edges (interactions) in the network.
        initial_nodes (list): A list storing the initial list of nodes.
        resources (pd.DataFrame): A DataFrame storing the interactions from the Omnipath database.
        ontology (Ontology): An instance of the Ontology class.

    Methods:
        copy(): Returns a deep copy of the Network instance.
        add_node(node: str): Adds a node to the network.
        remove_node(node: str): Removes a node from the network.
        add_edge(edge: pd.DataFrame): Adds an edge to the network.
        load_network_from_sif(sif_file): Loads a network from a SIF (Simple Interaction Format) file.
        add_paths_to_edge_list(paths): Adds paths to the edge list of the network.
        connect_nodes(only_signed: bool = False, consensus_only: bool = False): Connects all the nodes in the network.
        is_connected(): Checks if all the nodes in the network are connected.
        filter_unsigned_paths(paths: list[tuple], consensus: bool): Filters out unsigned paths from the provided list of paths.
        connect_subgroup(group: (str | pd.DataFrame | list[str]), maxlen: int = 1, only_signed: bool = False, consensus: bool = False): Connects all the nodes in a particular subgroup.
        complete_connection(maxlen: int = 2, minimal: bool = True, k_mean: Literal['tight', 'extensive'] = 'tight', only_signed: bool = False, consensus: bool = False, connect_node_when_first_introduced: bool = True): Connects all nodes of a network object.
        connect_component(comp_A: (str | pd.DataFrame | list[str]), comp_B: (str | pd.DataFrame | list[str]), maxlen: int = 2, mode: str = 'OUT', only_signed: bool = False, compress: bool = False): Connects subcomponents of a network object.
        convert_edgelist_into_genesymbol(): Converts the edge dataframe from uniprot to genesymbol.
        connect_genes_to_phenotype(phenotype: str = None, id_accession: str = None, sub_genes: list[str] = None, maxlen: int = 2, only_signed: bool = False, compress: bool = False): Connects genes to a phenotype.
    """

    def __init__(self,
                 initial_nodes: list[str] = None,
                 sif_file=None,
                 resources=None):
        self.nodes = pd.DataFrame(columns=["Genesymbol", "Uniprot", "Type"])
        self.edges = pd.DataFrame(columns=["source", "target", "Type", "Effect", "References"])
        self.initial_nodes = initial_nodes
        self.ontology = Ontology()
        if resources is not None and isinstance(resources, pd.DataFrame) and not resources.empty:
            self.resources = resources.copy()
        else:
            res = Resources()
            res.load_all_omnipath_interactions()
            self.resources = res.interactions.copy()
        if initial_nodes:
            for node in initial_nodes:
                self.add_node(node)
            self.drop_missing_nodes()
            self.nodes.reset_index(inplace=True, drop=True)
        elif sif_file:
            self.initial_nodes = []
            self.load_network_from_sif(sif_file)

        self.connect = Connections(self.resources)

    def copy(self):
        new_instance = copy.deepcopy(self)
        return new_instance

    def check_nodes(self, nodes: list[str]) -> list[str]:
        """
        This function checks if the nodes exist in the resources database and returns the nodes that are present.

        Parameters:
        - nodes: A list of node identifiers (strings). These are the nodes to be checked.

        Returns:
        - A list of node identifiers (strings) that exist in the resources database. If a node from the input list does not exist in the resources database, it is not included in the output list.

        The function works by iterating over the input list of nodes. For each node, it checks if the node exists in the 'source' or 'target' columns of the resources database. If the node exists, it is added to the output list.
        """
        return [node for node in nodes if
                node in self.resources["source"].unique() or node in self.resources["target"].unique()]

    def drop_missing_nodes(self):
        """
        This function drops the nodes that are not present in the resources database and print a warning with the name of the missing nodes.

        The function works as follows:
        1. It first calls the `check_nodes` function to get a list of nodes that exist in the resources' database.
        2. It then finds the nodes in the network that are not in this list, and removes them from the network.
        3. If there are any missing nodes, it prints a warning with their names.

        This function does not return anything. It modifies the `nodes` attribute of the `Network` object in-place.
        """
        # Get the list of nodes that exist in the resources database
        existing_nodes = self.check_nodes(self.nodes["Uniprot"].tolist())

        # Find the nodes in the network that are not in the list of existing nodes
        missing_nodes = [node for node in self.nodes["Uniprot"].tolist() if node not in existing_nodes]

        # Remove the missing nodes from the network
        self.nodes = self.nodes[~self.nodes["Uniprot"].isin(missing_nodes)]

        # Print a warning with the name of the missing nodes
        if missing_nodes:
            print(
                "Warning: The following nodes were not found in the resources database and have been removed from the "
                "network:",
                ", ".join(missing_nodes))
        return

    def add_node(self, node: str, from_sif: bool = False):
        """
        Adds a node to the network. The node is added to the nodes DataFrame of the network. The function checks the syntax
        for the genesymbol to ensure it is correct. If the node is a complex, it is added with the 'Genesymbol' as the complex
        string and 'Uniprot' as the node. Otherwise, it is added with the 'Genesymbol' as the genesymbol and 'Uniprot' as the
        uniprot. The 'Type' is set as 'NaN' for all new nodes.

        Parameters:
        - node: A string representing the node to be added. The node can be represented by either its Genesymbol or Uniprot identifier.

        Returns:
        None. The function modifies the network object in-place by adding the node to the nodes DataFrame.
        """

        if from_sif:
            #check that the new entry node can be translated using the function mapping node identifier (all the output of the function should be None)
            #if it cannot be translated, print an error message but add the node to the network anyway

            complex_string, genesymbol, uniprot = mapping_node_identifier(node)
            if not complex_string and not genesymbol and not uniprot:
                print("Error: node %s could not be automatically translated" % node)
                new_entry = {"Genesymbol": node, "Uniprot": node, "Type": "NaN"}
                self.nodes.loc[len(self.nodes)] = new_entry
                self.nodes = self.nodes.drop_duplicates().reset_index(drop=True)

            new_entry = {"Genesymbol": genesymbol, "Uniprot": uniprot, "Type": "NaN"}
            self.nodes.loc[len(self.nodes)] = new_entry
            self.nodes = self.nodes.drop_duplicates().reset_index(drop=True)
            self.initial_nodes.append(new_entry["Genesymbol"])
            self.initial_nodes = list(set(self.initial_nodes))
            return

        complex_string, genesymbol, uniprot = mapping_node_identifier(node)

        if complex_string:
            new_entry = {"Genesymbol": complex_string, "Uniprot": node, "Type": "NaN"}
        else:
            new_entry = {"Genesymbol": genesymbol, "Uniprot": uniprot, "Type": "NaN"}

        if not self.check_node_existence(uniprot):
            print("Error: node %s is not present in the resources database" % node)
            return

        self.nodes.loc[len(self.nodes)] = new_entry

        self.nodes = self.nodes.drop_duplicates().reset_index(drop=True)
        return

    def remove_node(self, node: str):
        """
        Removes a node from the network. The node is removed from both the list of nodes and the list of edges.

        Parameters:
        - node: A string representing the node to be removed. The node can be represented by either its Genesymbol or Uniprot identifier.

        Returns:
        None. The function modifies the network object in-place by removing the node from the nodes DataFrame and any associated edges from the edges DataFrame.
        """
        # Remove the node from the nodes DataFrame
        self.nodes = self.nodes[(self.nodes.Genesymbol != node) & (self.nodes.Uniprot != node)]

        # Translate the node identifier to Uniprot
        node = mapping_node_identifier(node)[2]

        # Remove any edges associated with the node from the edges DataFrame
        self.edges = self.edges[~self.edges[['source', 'target']].isin([node]).any(axis=1)]

        return

    def add_edge(self, edge: pd.DataFrame):
        """
        This method adds an interaction to the list of interactions while converting it to the Omniflow-network format.
        It checks if the edge represents inhibition or stimulation and sets the effect accordingly. It also checks if the
        nodes involved in the interaction are already present in the network, if not, it adds them.

        Parameters:
        - edge: A pandas DataFrame representing the interaction. The DataFrame should contain columns for 'source', 'target',
                 'type', and 'references'. The 'source' and 'target' columns represent the nodes involved in the interaction.
                 The 'type' column represents the type of interaction. The 'references' column contains the references for the interaction.

        Returns:
        None. The function modifies the network object in-place by adding the interaction to the edges DataFrame and adding
        any new nodes to the nodes DataFrame.
        """

        # Check if the edge represents inhibition or stimulation and set the effect accordingly
        effect = check_sign(edge)
        references = edge["references"].values[0]
        # Get the type value from the edge DataFrame or set it to None
        edge_type = edge["type"].values[0] if "type" in edge.columns else None

        # Create a new DataFrame with edge information, including handling None for type
        df_edge = pd.DataFrame({
            "source": edge["source"],
            "target": edge["target"],
            "Type": edge_type,
            "Effect": effect,
            "References": references
        })

        # Convert the "Uniprot" column to a set for efficient membership test
        uniprot_nodes = set(self.nodes["Uniprot"].unique())

        # add the new nodes to the nodes dataframe
        if edge["source"].values[0] not in uniprot_nodes:
            self.add_node(edge["source"].values[0])
        if edge["target"].values[0] not in uniprot_nodes:
            self.add_node(edge["target"].values[0])

        # if in the edge dataframe there is an edge with the same source, target and effect, merge the references
        existing_edge = self.edges[(self.edges["source"] == edge["source"].values[0]) &
                                   (self.edges["target"] == edge["target"].values[0]) &
                                   (self.edges["Effect"] == effect)]
        if not existing_edge.empty:
            self.edges.loc[existing_edge.index, "References"] += "; " + str(references)
        else:
            # Concatenate the new edge DataFrame with the existing edges in the graph
            self.edges = pd.concat([self.edges, df_edge])

        self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def remove_edge(self, node1: str, node2: str):
        """
        This function removes an edge from the network. It takes the source node and target node as input and removes
        the edge from the edges DataFrame.

        Parameters:
        - node1: A string representing the source node of the edge.
        - node2: A string representing the target node of the edge.

        Returns:
        None. The function modifies the network object in-place by removing the edge from the edges DataFrame.
        """
        # check if node1 and node2 are in genesymbol format or uniprot format
        if check_gene_list_format([node1]):
            node1 = mapping_node_identifier(node1)[2]
        if check_gene_list_format([node2]):
            node2 = mapping_node_identifier(node2)[2]

        # Remove the edge from the edges DataFrame, if the effect or the nodes are not present, print a warning
        if not self.edges[(self.edges["source"] == node1) & (self.edges["target"] == node2)].empty:
            self.edges = self.edges[~((self.edges["source"] == node1) & (self.edges["target"] == node2))]
        else:
            print("Warning: The edge does not exist in the network, check syntax for ",
                  mapping_node_identifier(node1)[1], " and ", mapping_node_identifier(node2)[1])
        return

    def remove_disconnected_nodes(self):
        """
        This function removes nodes from the network that are not connected to any other nodes.

        Returns:
        None. The function modifies the network object in-place by removing the disconnected nodes from the nodes DataFrame.
        """
        # Get the list of nodes that are not connected to any other nodes
        disconnected_nodes = self.nodes[~self.nodes["Uniprot"].isin(self.edges["source"]) &
                                        ~self.nodes["Uniprot"].isin(self.edges["target"])]

        # Remove the disconnected nodes from the nodes DataFrame
        self.nodes = self.nodes[~self.nodes["Uniprot"].isin(disconnected_nodes["Uniprot"])]

        return

    def modify_node_name(self, old_name: str, new_name: str,
                         type: Literal['Genesymbol', 'Uniprot', 'both'] = 'Genesymbol'):
        """
        This function modifies the name of a node in the network. It takes the old name of the node and the new name
        as input and modifies the name of the node in the nodes and in the edges DataFrame. If type is set to 'Genesymbol',
        it modifies the genesymbol name of the node in the nodes DataFrame. If type is set to 'Uniprot', it modifies the
         uniprot name of the node in the edges DataFrame. If type is set to 'both', it modifies both the genesymbol and
            uniprot names of the node in the nodes and edges DataFrame.


        Parameters:
        - old_name: A string representing the old name of the node.
        - new_name: A string representing the new name of the node.
        - type: A string indicating the type of name to be modified. It can be 'Genesymbol', 'Uniprot', or 'both'. Default is 'Genesymbol'.

        Returns:
        None. The function modifies the network object in-place by changing the name of the node in the nodes DataFrame.
        """

        if type == 'Genesymbol':
            self.nodes.loc[self.nodes["Genesymbol"] == old_name, "Genesymbol"] = new_name
        elif type == 'Uniprot':
            self.nodes.loc[self.nodes["Uniprot"] == old_name, "Uniprot"] = new_name
            # Update the source and target columns in the edges DataFrame
            self.edges.loc[self.edges["source"] == old_name, "source"] = new_name
            self.edges.loc[self.edges["target"] == old_name, "target"] = new_name
        elif type == 'both':
            self.nodes.loc[self.nodes["Genesymbol"] == old_name, "Genesymbol"] = new_name
            # check if it is possible to translate the genesymbol to uniprot
            try:
                new_name_uniprot = mapping_node_identifier(new_name)[2]
                old_name_uniprot = mapping_node_identifier(old_name)[2]
            except:
                new_name_uniprot = new_name
                old_name_uniprot = old_name
            self.nodes.loc[self.nodes["Uniprot"] == old_name_uniprot, "Uniprot"] = new_name_uniprot
            # Update the source and target columns in the edges DataFrame
            self.edges.loc[self.edges["source"] == old_name_uniprot, "source"] = new_name_uniprot
            self.edges.loc[self.edges["target"] == old_name_uniprot, "target"] = new_name_uniprot
        else:
            print("Error: Invalid type. Please choose 'Genesymbol', 'Uniprot', or 'both'.")

        return

    def print_my_paths(self, node1: str, node2: str, maxlen: int = 2, genesymbol: bool = True):
        """
        This function prints all the paths between two nodes in the network. It uses the `find_paths` method from the
        `Connections` class to find all the paths between the two nodes in the Network object. If no paths are found,
        it prints a warning message. If one of the selected nodes is not present in the network, it prints an error message.

        Parameters:
        - node1: A string representing the source node.
        - node2: A string representing the target node.
        - maxlen: An integer representing the maximum length of the paths to be searched for. Default is 2.
        - genesymbol: A boolean flag indicating whether to print the paths in genesymbol format. Default is True.

        Returns:
        None. The function prints the paths between the two nodes in the network.
        """

        # check if node1 and node2 are in genesymbol format or uniprot format
        if check_gene_list_format([node1]):
            node1 = mapping_node_identifier(node1)[2]
        if check_gene_list_format([node2]):
            node2 = mapping_node_identifier(node2)[2]

        # Check if the nodes exist in the network
        if node1 not in self.nodes["Uniprot"].tolist() or node2 not in self.nodes["Uniprot"].tolist():
            print("Error: One or both of the selected nodes are not present in the network.")
            return
        connect = Connections(self.edges)
        paths = connect.find_paths(node1, node2, maxlen=maxlen)

        if not paths:
            print("Warning: No paths found between source: ", node1, " and target: ", node2)
            return

        if genesymbol:
            paths = translate_paths(paths)

        # Print all the paths
        for path in paths:
            print(path)

        return

    def remove_path(self, path: list[str]):
        """
        This function removes a path from the network. It takes a list of nodes representing the path and removes all the edges
        between the nodes in the path.

        Parameters:
        - path: A list of nodes representing the path to be removed. The nodes can be represented as strings or tuples.

        Returns:
        None. The function modifies the network object in-place by removing the edges between the nodes in the path.
        """
        # check if node1 and node2 are in genesymbol format or uniprot format
        if check_gene_list_format(path):
            path = [mapping_node_identifier(node)[2] for node in path]

        # Iterate through the nodes in the path and remove the edges between them
        for i in range(0, len(path)):
            if i == len(path) - 1:
                break
            self.remove_edge(path[i], path[i + 1])
        return

    def load_network_from_sif(self, sif_file):
        """
        Load a network object from a SIF (Simple Interaction Format) file.

        Parameters:
        - sif_file: A string representing the path to the SIF file.

        Returns:
        None. The function modifies the network object in-place.
        """
        interactions = []
        node_set = set()

        def determine_effect(interaction_type):
            """
            Determine the effect based on the interaction type.

            Parameters:
            - interaction_type: A string representing the type of interaction.

            Returns:
            - A string representing the effect of the interaction. If the interaction type is not recognized, it returns "undefined".
            """
            effect_types = {
                "1": "stimulation",
                "activate": "stimulation",
                "stimulate": "stimulation",
                "phosphorilate": "undefined",
                "stimulation": "stimulation",
                "->": "stimulation",
                "-|": "inhibition",
                "-1": "inhibition",
                "inhibit": "inhibition",
                "block": "inhibition",
                "inhibition": "inhibition",
                "form_complex": "form complex",
                "form-complex": "form complex",
                "complex_formation": "form complex"
            }
            return effect_types.get(interaction_type, "undefined")

        with open(sif_file, "r") as f:
            for line in f:
                if line.startswith('#'):  # Skip comment lines
                    continue
                interaction = line.strip().split()
                if len(interaction) < 3:
                    continue  # Skip malformed lines

                effect = determine_effect(interaction[1])

                interactions.append({
                    "source": mapping_node_identifier(interaction[0])[2] if check_gene_list_format(
                        [interaction[0]]) else interaction[0],
                    "target": mapping_node_identifier(interaction[2])[2] if check_gene_list_format(
                        [interaction[2]]) else interaction[2],
                    "Type": interaction[3] if len(interaction) > 3 else None,
                    "Effect": effect,
                    "References": "SIF file"
                })
                node_set.update([interaction[0], interaction[2]])

        # Create or update the edges DataFrame
        df_edge = pd.DataFrame(interactions)
        self.edges = pd.concat([self.edges, df_edge], ignore_index=True)

        # Update the nodes list
        nodes = list(node_set)
        for node in nodes:
            self.add_node(node, from_sif=True)

        return

    def add_paths_to_edge_list(self, paths):
        """
        This method adds paths to the edge list of the network. A path is a sequence of nodes where each node is
        connected to the next node in the sequence. The function checks if there is an interaction between each pair of nodes
        in the path in the resources' database. If an interaction exists, it is added to the edge list of the network.

        Parameters:
        - paths: A list of paths, where each path is a sequence of nodes. A node can be a string or a tuple.

        Returns:
        None. The function modifies the network object in-place by adding the interactions to the edges DataFrame.
        """
        # Access the resources database
        database = self.resources

        # Iterate through the list of paths
        for path in paths:
            # Handle single string or tuple
            if isinstance(path, (str)):
                path = [path]

            # Iterate through the nodes in the path
            for i in range(0, len(path)):
                # Break the loop if it's the last node in the path
                if i == len(path) - 1:
                    break

                # Check if there is an interaction between the current node and the next node in the resources database
                interaction = database.loc[(database["source"] == path[i]) &
                                           (database["target"] == path[i + 1])]

                # If an interaction exists, add it to the edge list of the network
                if not interaction.empty:
                    if not ((self.edges['source'] == interaction['source'].values[0]) & (
                        self.edges['target'] == interaction['target'].values[0])).any():
                        self.add_edge(interaction)

        # Remove duplicate edges from the edge list
        self.edges = self.edges.drop_duplicates().reset_index(drop=True)

        return

    def add_cascade_to_edge_list(self, cascades):
        """
        This function adds cascades to the edge list of the network. A cascade is a sequence of nodes where each node is
        connected to the next node in the sequence. The function checks if there is an interaction between each pair of nodes
        in the cascade in the resources' database. If an interaction exists, it is added to the edge list of the network.

        Parameters:
        - cascades: A list of cascades, where each cascade is a sequence of nodes.

        Returns:
        None. The function modifies the network object in-place.
        """
        database = self.resources

        for cascade in cascades:
            interaction_in = database.loc[(database["source"] == cascade[0]) &
                                          (database["target"] == cascade[1])]
            if interaction_in.empty:
                print("Empty interaction for node ", cascade[0], " and ", cascade[1])
            else:
                self.add_edge(interaction_in)
        self.edges = self.edges.drop_duplicates().reset_index(drop=True)

        return

    def connect_nodes(self,
                      only_signed: bool = False,
                      consensus_only: bool = False):
        """
        Basic node connections. It adds all the interactions found in the omnipath database.
        Once an interaction is found it will be added to the list of edges.
        The only_signed flag makes sure that just signed interaction will be added to the network, while "consensus_only"
        makes sure that just signed interaction with consensus among references will be included.

        Parameters:
        - only_signed: A boolean flag indicating whether to only add signed interactions to the network.
        - consensus_only: A boolean flag indicating whether to only add signed interactions with consensus among references to the network.

        Returns:
        None. The function modifies the network object in-place.
        """

        if len(self.nodes) == 1:
            print("Number of node insufficient to create connection")
            return

        def add_edge_if_not_empty_and_signed(node1, node2):
            """
            Helper function to add an edge to the network if the interaction is not empty and, if the `only_signed` flag is set, the interaction is signed.

            Parameters:
            - node1: The source node of the interaction.
            - node2: The target node of the interaction.

            Returns:
            None. The function modifies the network object in-place.
            """
            interaction = self.resources.loc[(self.resources["source"] == node1) &
                                             (self.resources["target"] == node2)]
            if not interaction.empty and (not only_signed or check_sign(interaction, consensus_only) != "undefined"):
                self.add_edge(interaction)

        for node1, node2 in combinations(self.nodes["Uniprot"], 2):
            add_edge_if_not_empty_and_signed(node1, node2)
            add_edge_if_not_empty_and_signed(node2, node1)
        return

    def is_connected(self) -> bool:
        """
        Return True if all the nodes in the nodes list in the Network object are connected, otherwise it returns False
        """
        # Create a set to store visited nodes
        visited = set()

        # Function for Depth-First Search
        def dfs(node):
            visited.add(node)
            for _, edge in self.edges.iterrows():
                if edge['source'] == node and edge['target'] not in visited:
                    dfs(edge['target'])
                elif edge['target'] == node and edge['source'] not in visited:
                    dfs(edge['source'])

        # Start DFS from the first node
        dfs(self.nodes.iloc[0]['Uniprot'])

        # Check if all nodes are visited
        return set(self.nodes['Uniprot']) == visited

    def filter_unsigned_paths(self, paths: list[tuple], consensus: bool) -> list[tuple]:
        """
        This function filters out unsigned paths from the provided list of paths. An unsigned path is a path where at least
        one interaction does not have a defined sign (stimulation or inhibition). The function checks each interaction in each
        path, and if the interaction is unsigned, the path is not included in the output list.

        Parameters:
        - paths: A list of paths, where each path is a sequence of nodes.
        - consensus: A boolean indicating whether to check for consensus among references when determining the sign of an interaction.

        Returns:
        - A list of paths where all interactions in each path are signed.
        """

        interactions = self.resources
        filtered_paths = []
        for path in paths:
            if isinstance(path, (str)):
                path = [path]
            is_full_signed = True
            for i in range(0, len(path)):
                if i == len(path) - 1:
                    break
                interaction = interactions.loc[(interactions["source"] == path[i]) &
                                               (interactions["target"] == path[i + 1])]
                if not interaction.empty and check_sign(interaction, consensus) == "undefined":
                    is_full_signed = False
                    break
            if is_full_signed:
                filtered_paths.append(path)

        return filtered_paths

    def check_node_existence(self, node: str) -> bool:
        """
        This function checks if a node exists in the resources' database.

        Parameters:
        - node: A string representing the node to be checked.

        Returns:
        - A boolean indicating whether the node exists in the resources' database.
        """
        return node in self.resources["source"].unique() or node in self.resources["target"].unique()

    def connect_subgroup(self,
                         group: (str | pd.DataFrame | list[str]),
                         maxlen: int = 1,
                         only_signed: bool = False,
                         consensus: bool = False
                         ):
        """
        This function is used to connect all the nodes in a particular subgroup. It iterates over all pairs of nodes in the
        subgroup and finds paths between them in the resources' database. If a path is found, it is added to the edge list of
        the network. The function also filters out unsigned paths if the `only_signed` flag is set to True.

        Parameters:
        - group: A list of nodes representing the subgroup to connect. Nodes can be represented as strings, pandas DataFrame, or list of strings.
        - maxlen: The maximum length of the paths to be searched for in the resources' database. Default is 1.
        - only_signed: A boolean flag indicating whether to only add signed interactions to the network. Default is False.
        - consensus: A boolean flag indicating whether to only add signed interactions with consensus among references to the network. Default is False.

        Returns:
        None. The function modifies the network object in-place.
        """

        if not check_gene_list_format(group):
            uniprot_gene_list = group
        else:
            uniprot_gene_list = [mapping_node_identifier(i)[2] for i in group]
        if len(uniprot_gene_list) == 1:
            print("Number of node insufficient to create connection")
        else:
            for node1, node2 in combinations(uniprot_gene_list, 2):
                i = 0
                paths_in = []
                paths_out = []
                while i <= maxlen:
                    if not paths_out:
                        paths_out = self.connect.find_paths(node1, node2, maxlen=i)
                        if only_signed:
                            paths_out = self.filter_unsigned_paths(paths_out, consensus)
                    if not paths_in:
                        paths_in = self.connect.find_paths(node2, node1, maxlen=i)
                        if only_signed:
                            paths_in = self.filter_unsigned_paths(paths_in, consensus)
                    if not paths_in or not paths_out and i <= maxlen:
                        i += 1
                    if (paths_in or paths_out) and i > maxlen or (paths_in and paths_out):
                        paths = paths_out + paths_in
                        self.add_paths_to_edge_list(paths)
                        break
        return

    def complete_connectionv2(self,
                              maxlen: int = 2,
                              minimal: bool = True,
                              k_mean: Literal['tight', 'extensive'] = 'tight',
                              only_signed: bool = False,
                              consensus: bool = False,
                              connect_node_when_first_introduced: bool = True):
        """
        This function attempts to connect all nodes of a network object using one of the methods presented in the Connection
        object. This is a core characteristic of this package and the user should have the possibility to choose different
        methods to enrich its Network object.

        Parameters:
        - maxlen: The maximum length of the paths to be searched for. Default is 2.
        - minimal: A boolean flag indicating whether to reset the object connect_network, updating the possible list of
          paths. Default is True.
        - k_mean: The search mode, which can be 'tight' or 'extensive'. Default is 'tight'.
        - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
        - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
        - connect_node_when_first_introduced: A boolean flag indicating whether to connect nodes when first introduced.
          Default is True.

        Returns:
        None. The function modifies the network object in-place.
        """

        # Set the search depth based on the k_mean parameter
        if k_mean == 'tight':
            i_search = 3
        elif k_mean == 'extensive':
            i_search = 4

        # Copy the nodes
        nodes = self.nodes.copy()

        # Create a Connections object for the edges
        connect_network = Connections(self.edges)

        # Iterate through all combinations of nodes
        for node1, node2 in combinations(nodes["Uniprot"], 2):
            if not self.check_node_existence(node1) or not self.check_node_existence(node2):
                print(
                    "Error: node %s is not present in the resources database" % node1 if not self.check_node_existence(
                        node1) else "Error: node %s is not present in the resources database" % node2)
                continue
            i = 0
            # Reset the object connect_network, updating the possible list of paths if minimal is True
            if minimal:
                connect_network = Connections(self.edges)

            # As first step, make sure that there is at least one path between two nodes in the network separated by len = maxlen
            paths_in = connect_network.find_shortest_path(node2, node1, max_len=maxlen)
            paths_out = connect_network.find_shortest_path(node1, node2, max_len=maxlen)

            if not paths_in:
                paths_in = self.connect.find_all_shortest_paths(node2, node1, max_len=i_search)
            if only_signed and paths_in:
                paths_in = self.filter_unsigned_paths(paths_in, consensus)
            if paths_in:
                self.add_paths_to_edge_list(paths_in)
            if connect_node_when_first_introduced:
                self.connect_nodes(only_signed, consensus)
                self.edges = self.edges.drop_duplicates().reset_index(drop=True)

            if not paths_out:
                paths_out = self.connect.find_all_shortest_paths(node1, node2, max_len=i_search)
            if only_signed and paths_out:
                paths_out = self.filter_unsigned_paths(paths_out, consensus)
            if paths_out:
                self.add_paths_to_edge_list(paths_out)
            if connect_node_when_first_introduced:
                self.connect_nodes(only_signed, consensus)
                self.edges = self.edges.drop_duplicates().reset_index(drop=True)

        # If connect_node_when_first_introduced is False, connect nodes after all paths have been found
        if not connect_node_when_first_introduced:
            self.connect_nodes(only_signed, consensus)
            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def complete_connectionv3(self,
                              maxlen: int = 2,
                              minimal: bool = True,
                              only_signed: bool = False,
                              consensus: bool = False,
                              connect_node_when_first_introduced: bool = True):
        """
        This function attempts to connect all nodes of a network object using one of the methods presented in the Connection
        object. This is a core characteristic of this package and the user should have the possibility to choose different
        methods to enrich its Network object.

        Parameters:
        - maxlen: The maximum length of the paths to be searched for. Default is 2.
        - minimal: A boolean flag indicating whether to reset the object connect_network, updating the possible list of
          paths. Default is True.
        - k_mean: The search mode, which can be 'tight' or 'extensive'. Default is 'tight'.
        - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
        - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
        - connect_node_when_first_introduced: A boolean flag indicating whether to connect nodes when first introduced.
          Default is True.

        Returns:
        None. The function modifies the network object in-place.
        """

        # Copy the nodes
        nodes = self.nodes.copy()

        # Create a Connections object for the edges
        connect_network = Connections(self.edges)

        # Iterate through all combinations of nodes
        for node1, node2 in combinations(nodes["Uniprot"], 2):
            if not self.check_node_existence(node1) or not self.check_node_existence(node2):
                print(
                    "Error: node %s is not present in the resources database" % node1 if not self.check_node_existence(
                        node1) else "Error: node %s is not present in the resources database" % node2)
                continue
            i = 0
            # Reset the object connect_network, updating the possible list of paths if minimal is True
            if minimal:
                connect_network = Connections(self.edges)

            # As first step, make sure that there is at least one path between two nodes in the network separated by len = maxlen
            paths_in = connect_network.bfs(node2, node1)
            paths_out = connect_network.bfs(node1, node2)

            if not paths_in:
                i = 0
                while i <= maxlen:
                    paths_in = self.connect.find_paths(node2, node1, maxlen=i)
                    if only_signed:
                        paths_in = self.filter_unsigned_paths(paths_in, consensus)
                    if paths_in:
                        self.add_paths_to_edge_list(paths_in)
                        if connect_node_when_first_introduced:
                            self.connect_nodes(only_signed, consensus)
                            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
                        break
                    else:
                        i += 1

            if not paths_out:
                i = 0
                while i <= maxlen:
                    paths_out = self.connect.find_paths(node1, node2, maxlen=i)
                    if only_signed:
                        paths_out = self.filter_unsigned_paths(paths_out, consensus)
                    if paths_out:
                        self.add_paths_to_edge_list(paths_out)
                        if connect_node_when_first_introduced:
                            self.connect_nodes(only_signed, consensus)
                            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
                        break
                    else:
                        i += 1

        # If connect_node_when_first_introduced is False, connect nodes after all paths have been found
        if not connect_node_when_first_introduced:
            self.connect_nodes(only_signed, consensus)
            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def complete_connectionv5(self,
                              maxlen: int = 2,
                              minimal: bool = True,
                              only_signed: bool = False,
                              consensus: bool = False,
                              connect_node_when_first_introduced: bool = True):
        """
        This function attempts to connect all nodes of a network object using one of the methods presented in the Connection
        object. This is a core characteristic of this package and the user should have the possibility to choose different
        methods to enrich its Network object.

        Parameters:
        - maxlen: The maximum length of the paths to be searched for. Default is 2.
        - minimal: A boolean flag indicating whether to reset the object connect_network, updating the possible list of
          paths. Default is True.
        - k_mean: The search mode, which can be 'tight' or 'extensive'. Default is 'tight'.
        - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
        - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
        - connect_node_when_first_introduced: A boolean flag indicating whether to connect nodes when first introduced.
          Default is True.

        Returns:
        None. The function modifies the network object in-place.
        """

        # Copy the nodes
        nodes = self.nodes.copy()

        # Create a Connections object for the edges
        connect_network = Connections(self.edges)

        # Iterate through all combinations of nodes
        for node1, node2 in combinations(nodes["Uniprot"], 2):
            if not self.check_node_existence(node1) or not self.check_node_existence(node2):
                print(
                    "Error: node %s is not present in the resources database" % node1 if not self.check_node_existence(
                        node1) else "Error: node %s is not present in the resources database" % node2)
                continue
            i = 0
            # Reset the object connect_network, updating the possible list of paths if minimal is True
            if minimal:
                connect_network = Connections(self.edges)

            # As first step, make sure that there is at least one path between two nodes in the network separated by len = maxlen
            paths_in = connect_network.find_shortest_path(node2, node1, max_len=maxlen)
            paths_out = connect_network.find_shortest_path(node1, node2, max_len=maxlen)

            if not paths_in:
                i = 0
                while i <= maxlen:
                    paths_in = self.connect.find_paths(node2, node1, maxlen=i)
                    if only_signed:
                        paths_in = self.filter_unsigned_paths(paths_in, consensus)
                    if paths_in:
                        self.add_paths_to_edge_list(paths_in)
                        if connect_node_when_first_introduced:
                            self.connect_nodes(only_signed, consensus)
                            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
                        break
                    else:
                        i += 1

            if not paths_out:
                i = 0
                while i <= maxlen:
                    paths_out = self.connect.find_paths(node1, node2, maxlen=i)
                    if only_signed:
                        paths_out = self.filter_unsigned_paths(paths_out, consensus)
                    if paths_out:
                        self.add_paths_to_edge_list(paths_out)
                        if connect_node_when_first_introduced:
                            self.connect_nodes(only_signed, consensus)
                            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
                        break
                    else:
                        i += 1

        # If connect_node_when_first_introduced is False, connect nodes after all paths have been found
        if not connect_node_when_first_introduced:
            self.connect_nodes(only_signed, consensus)
            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def complete_connectionv6(self,
                              minimal: bool = True,
                              only_signed: bool = False,
                              consensus: bool = False,
                              connect_node_when_first_introduced: bool = True):
        """
        This function attempts to connect all nodes of a network object using one of the methods presented in the Connection
        object. This is a core characteristic of this package and the user should have the possibility to choose different
        methods to enrich its Network object.

        Parameters:
        - maxlen: The maximum length of the paths to be searched for. Default is 2.
        - minimal: A boolean flag indicating whether to reset the object connect_network, updating the possible list of
          paths. Default is True.
        - k_mean: The search mode, which can be 'tight' or 'extensive'. Default is 'tight'.
        - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
        - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
        - connect_node_when_first_introduced: A boolean flag indicating whether to connect nodes when first introduced.
          Default is True.

        Returns:
        None. The function modifies the network object in-place.
        """

        # Copy the nodes
        nodes = self.nodes.copy()

        # Create a Connections object for the edges
        connect_network = Connections(self.edges)

        # Iterate through all combinations of nodes
        for node1, node2 in combinations(nodes["Uniprot"], 2):
            if not self.check_node_existence(node1) or not self.check_node_existence(node2):
                print(
                    "Error: node %s is not present in the resources database" % node1 if not self.check_node_existence(
                        node1) else "Error: node %s is not present in the resources database" % node2)
                continue

            # Reset the object connect_network, updating the possible list of paths if minimal is True
            if minimal:
                connect_network = Connections(self.edges)

            # As first step, make sure that there is at least one path between two nodes in the network separated by len = maxlen
            paths_in = connect_network.bfs(node2, node1)
            paths_out = connect_network.bfs(node1, node2)

            if not paths_in:
                paths_in = self.connect.bfs(node2, node1)
                if only_signed:
                    paths_in = self.filter_unsigned_paths(paths_in, consensus)
                if paths_in:
                    self.add_paths_to_edge_list(paths_in)
                    if connect_node_when_first_introduced:
                        self.connect_nodes(only_signed, consensus)
                        self.edges = self.edges.drop_duplicates().reset_index(drop=True)

            if not paths_out:
                paths_out = self.connect.bfs(node1, node2)
                if only_signed:
                    paths_out = self.filter_unsigned_paths(paths_out, consensus)
                if paths_out:
                    self.add_paths_to_edge_list(paths_out)
                    if connect_node_when_first_introduced:
                        self.connect_nodes(only_signed, consensus)
                        self.edges = self.edges.drop_duplicates().reset_index(drop=True)

        # If connect_node_when_first_introduced is False, connect nodes after all paths have been found
        if not connect_node_when_first_introduced:
            self.connect_nodes(only_signed, consensus)
            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def complete_connectionv4(self,
                              maxlen: int = 2,
                              minimal: bool = True,
                              only_signed: bool = False,
                              consensus: bool = False,
                              connect_node_when_first_introduced: bool = True):
        """
        This function attempts to connect all nodes of a network object using one of the methods presented in the Connection
        object. This is a core characteristic of this package and the user should have the possibility to choose different
        methods to enrich its Network object.

        Parameters:
        - maxlen: The maximum length of the paths to be searched for. Default is 2.
        - minimal: A boolean flag indicating whether to reset the object connect_network, updating the possible list of
          paths. Default is True.
        - k_mean: The search mode, which can be 'tight' or 'extensive'. Default is 'tight'.
        - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
        - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
        - connect_node_when_first_introduced: A boolean flag indicating whether to connect nodes when first introduced.
          Default is True.

        Returns:
        None. The function modifies the network object in-place.
        """

        # Copy the nodes
        nodes = self.nodes.copy()

        # Create a Connections object for the edges
        connect_network = Connections(self.edges)

        # Iterate through all combinations of nodes
        for node1, node2 in combinations(nodes["Uniprot"], 2):
            if not self.check_node_existence(node1) or not self.check_node_existence(node2):
                print(
                    "Error: node %s is not present in the resources database" % node1 if not self.check_node_existence(
                        node1) else "Error: node %s is not present in the resources database" % node2)
                continue
            i = 0
            # Reset the object connect_network, updating the possible list of paths if minimal is True
            if minimal:
                connect_network = Connections(self.edges)

            # As first step, make sure that there is at least one path between two nodes in the network separated by len = maxlen
            paths_in = connect_network.bfs(node2, node1)
            paths_out = connect_network.bfs(node1, node2)

            if not paths_in:
                paths_in = self.connect.find_all_shortest_paths(node2, node1, max_len=maxlen)
            if only_signed and paths_in:
                paths_in = self.filter_unsigned_paths(paths_in, consensus)
            if paths_in:
                self.add_paths_to_edge_list(paths_in)
                if connect_node_when_first_introduced:
                    self.connect_nodes(only_signed, consensus)
                    self.edges = self.edges.drop_duplicates().reset_index(drop=True)

            if not paths_out:
                paths_out = self.connect.find_all_shortest_paths(node1, node2, max_len=maxlen)
            if only_signed and paths_out:
                paths_out = self.filter_unsigned_paths(paths_out, consensus)
            if paths_out:
                self.add_paths_to_edge_list(paths_out)
                if connect_node_when_first_introduced:
                    self.connect_nodes(only_signed, consensus)
                    self.edges = self.edges.drop_duplicates().reset_index(drop=True)

        # If connect_node_when_first_introduced is False, connect nodes after all paths have been found
        if not connect_node_when_first_introduced:
            self.connect_nodes(only_signed, consensus)
            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def complete_connection(self,
                            maxlen: int = 2,
                            minimal: bool = True,
                            k_mean: Literal['tight', 'extensive'] = 'tight',
                            only_signed: bool = False,
                            consensus: bool = False,
                            connect_node_when_first_introduced: bool = True):
        """
        This function attempts to connect all nodes of a network object using one of the methods presented in the Connection
        object. This is a core characteristic of this package and the user should have the possibility to choose different
        methods to enrich its Network object.

        Parameters:
        - maxlen: The maximum length of the paths to be searched for. Default is 2.
        - minimal: A boolean flag indicating whether to reset the object connect_network, updating the possible list of
          paths. Default is True.
        - k_mean: The search mode, which can be 'tight' or 'extensive'. Default is 'tight'.
        - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
        - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
        - connect_node_when_first_introduced: A boolean flag indicating whether to connect nodes when first introduced.
          Default is True.

        Returns:
        None. The function modifies the network object in-place.
        """

        # Set the search depth based on the k_mean parameter
        if k_mean == 'tight':
            i_search = 3
        elif k_mean == 'extensive':
            i_search = 4

        # Copy the nodes
        nodes = self.nodes.copy()

        # Create a Connections object for the edges
        connect_network = Connections(self.edges)

        # Iterate through all combinations of nodes
        for node1, node2 in combinations(nodes["Uniprot"], 2):
            if not self.check_node_existence(node1) or not self.check_node_existence(node2):
                print(
                    "Error: node %s is not present in the resources database" % node1 if not self.check_node_existence(
                        node1) else "Error: node %s is not present in the resources database" % node2)
                continue
            i = 0
            # Reset the object connect_network, updating the possible list of paths if minimal is True
            if minimal:
                connect_network = Connections(self.edges)

            # Start a loop to find paths between nodes
            while i <= maxlen:
                # As first step, make sure that there is at least one path between two nodes in the network
                paths_in = connect_network.find_paths(node2, node1, maxlen=i)
                paths_out = connect_network.find_paths(node1, node2, maxlen=i)

                # If paths in both directions are found, break the loop
                if paths_in and paths_out:
                    break

                # If no paths are found and the maximum length has not been reached, increment the length and continue the loop
                if not (paths_in and paths_out) and i < maxlen:
                    i += 1
                    continue

                # If no paths are found and the maximum length has been reached, search for new paths in the database
                if not paths_in and i == maxlen:
                    flag = False
                    j = 0
                    # Start a loop to find new paths in the database
                    while not flag and j <= i_search:
                        paths_in = self.connect.find_paths(node2, node1, maxlen=j)
                        # Filter unsigned paths if only_signed is True
                        if only_signed:
                            paths_in = self.filter_unsigned_paths(paths_in, consensus)
                        # If no paths are found, increment the length and continue the loop
                        if not paths_in:
                            j += 1
                        else:
                            # If a path is found, add it to the edge list and connect nodes if connect_node_when_first_introduced is True
                            self.add_paths_to_edge_list(paths_in)
                            if connect_node_when_first_introduced:
                                self.connect_nodes(only_signed, consensus)
                                self.edges = self.edges.drop_duplicates().reset_index(drop=True)
                            flag = True

                # Repeat the same process for paths in the opposite direction
                if not paths_out and i == maxlen:
                    flag = False
                    j = 0
                    while not flag and j <= i_search:
                        paths_out = self.connect.find_paths(node1, node2, maxlen=j)
                        if only_signed:
                            paths_out = self.filter_unsigned_paths(paths_out, consensus)
                        if not paths_out:
                            j += 1
                        else:
                            self.add_paths_to_edge_list(paths_out)
                            if connect_node_when_first_introduced:
                                self.connect_nodes(only_signed, consensus)
                                self.edges = self.edges.drop_duplicates().reset_index(drop=True)
                            flag = True
                break

        # If connect_node_when_first_introduced is False, connect nodes after all paths have been found
        if not connect_node_when_first_introduced:
            self.connect_nodes(only_signed, consensus)
            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def connect_component(self,
                          comp_A: (str | list[str]),
                          comp_B: (str | list[str]),
                          maxlen: int = 2,
                          mode: Literal['OUT', 'IN', 'ALL'] = 'OUT',
                          only_signed: bool = False,
                          consensus: bool = False):
        """
        This function attempts to connect subcomponents of a network object using one of the methods presented in the Connection
        object. This is a core characteristic of this package and the user should have the possibility to choose different
        methods to enrich its Network object.

        Parameters:
        - comp_A: A string or list of strings representing the first component to connect.
        - comp_B: A string or list of strings representing the second component to connect.
        - maxlen: The maximum length of the paths to be searched for. Default is 2.
        - mode: The search mode, which can be 'OUT', 'IN', or 'ALL'. Default is 'OUT'.
        - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
        - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.

        Returns:
        None. The function modifies the network object in-place.

        """

        # Determine the search mode and find paths accordingly
        if mode == "IN":
            paths_in = self.connect.find_paths(comp_B, comp_A, maxlen=maxlen)
            paths = paths_in
        elif mode == "OUT":
            paths_out = self.connect.find_paths(comp_A, comp_B, maxlen=maxlen)
            paths = paths_out
        elif mode == "ALL":
            paths_out = self.connect.find_paths(comp_A, comp_B, maxlen=maxlen)
            paths_in = self.connect.find_paths(comp_B, comp_A, maxlen=maxlen)
            paths = paths_out + paths_in
        else:
            print("The only accepted modes are IN, OUT or ALL, please check the syntax")
            return

        # Filter unsigned paths if the only_signed flag is set
        if only_signed:
            paths = self.filter_unsigned_paths(paths, consensus)

        # Add the paths to the edge list
        self.add_paths_to_edge_list(paths)

        # Create sets of nodes for each component and the entire network
        all_nodes = set(self.nodes['Uniprot'].values)
        set_a = set(comp_A)
        set_b = set(comp_B)

        # Find the nodes that are not in either component
        set_c = all_nodes.difference(set_a)
        set_c = set_c.difference(set_b)
        set_c = list(set_c)

        # If there are nodes not in either component, connect them as a subgroup
        if len(set_c) > 0:
            self.connect_subgroup(set_c, only_signed=only_signed, maxlen=maxlen)

        return

    def connect_to_upstream_nodes(self,
                                  nodes_to_connect: Optional[List[str]] = None,
                                  depth: int = 1,
                                  rank: int = 1,
                                  only_signed: bool = True,
                                  consensus: bool = False) -> None:
        """
        Connects to upstream nodes based on the provided parameters.

        Parameters:
        - nodes_to_connect: A list of nodes to connect. If not provided, all nodes in the network are considered.
        - depth: The depth of the search for upstream nodes.
        - rank: The rank of the search for upstream nodes.

        Returns:
        None. The function modifies the network object in-place.
        """
        try:
            if nodes_to_connect is None:
                nodes_to_connect = self.nodes["Uniprot"].tolist()

            cascades = self.connect.find_upstream_cascades(nodes_to_connect, depth, rank)

            if only_signed:
                cascades = self.filter_unsigned_paths(cascades, consensus)
            self.add_cascade_to_edge_list(cascades)
            self.edges.drop_duplicates().reset_index(drop=True)
        except Exception as e:
            print(f"An error occurred while connecting to upstream nodes: {e}")
        return

    def convert_edgelist_into_genesymbol(self) -> pd.DataFrame:
        """
        This function converts the edge dataframe from uniprot to genesymbol.
        """

        def convert_identifier(x):
            identifiers = mapping_node_identifier(x)
            return identifiers[0] or identifiers[1]

        gs_edges = self.edges.copy()

        gs_edges["source"] = gs_edges["source"].apply(convert_identifier)
        gs_edges["target"] = gs_edges["target"].apply(convert_identifier)

        return gs_edges

    def connect_genes_to_phenotype(self,
                                   phenotype: str = None,
                                   id_accession: str = None,
                                   sub_genes: list[str] = None,
                                   maxlen: int = 2,
                                   only_signed: bool = False,
                                   compress: bool = False
                                   ):
        """
        This function connects genes to a phenotype based on the provided parameters. It retrieves phenotype markers,
        identifies unique Uniprot genes, and connects them to the network. It also has the option to compress the network
        by substituting specified genes with the phenotype name.

        Parameters:
        - phenotype: The phenotype to connect to. If not provided, it will be retrieved using the id_accession.
        - id_accession: The accession id of the phenotype. If not provided, the phenotype parameter must be given.
        - sub_genes: A list of genes to be considered for connection. If not provided, all nodes in the network are considered.
        - maxlen: The maximum length of the paths to be searched for.
        - only_signed: A boolean flag to indicate whether to filter unsigned paths.
        - compress: A boolean flag to indicate whether to substitute the specified genes with the phenotype name.

        Returns:
        None. The function modifies the network object in-place.
        """
        # Initialize lists for Uniprot and genesymbol genes
        uniprot_gene_list = []
        genesymbols_genes = []

        # Retrieve phenotype markers
        phenotype_genes = self.ontology.get_markers(phenotype=phenotype, id_accession=id_accession)
        if not phenotype_genes:
            print("Something went wrong while getting the markers for: ", phenotype, " and ", id_accession)
            print("Check URL and try again")
            return
        # Convert phenotype genes to Uniprot identifiers
        uniprot_genes = [mapping_node_identifier(i)[2] for i in phenotype_genes]
        # If sub_genes are provided, check their format and convert to Uniprot or genesymbol as needed
        if sub_genes:
            if check_gene_list_format(sub_genes):
                uniprot_gene_list = sub_genes
                genesymbols_genes = [mapping_node_identifier(i)[2] for i in sub_genes]
            else:
                uniprot_gene_list = [mapping_node_identifier(i)[2] for i in sub_genes]
                genesymbols_genes = sub_genes

        print("Starting connecting network's nodes to: ", phenotype_genes)
        # Identify unique Uniprot genes not already in the network
        unique_uniprot = set(uniprot_genes) - set(uniprot_gene_list if uniprot_gene_list else self.nodes["Uniprot"])
        # Identify unique genesymbols not already in the network
        unique_genesymbol = set(phenotype_genes) - set(
            genesymbols_genes if genesymbols_genes else self.nodes["Genesymbol"])
        # Connect the network's nodes to the unique Uniprot genes
        self.connect_component(uniprot_gene_list if uniprot_gene_list else self.nodes["Uniprot"].tolist(),
                               list(unique_uniprot),
                               mode="OUT",
                               maxlen=maxlen, only_signed=only_signed)

        # If compress is True, substitute specified genes with the phenotype name
        if compress:
            phenotype = phenotype or self.ontology.accession_to_phenotype_dict[id_accession]
            phenotype_modified = phenotype.replace(" ", "_")

            # Substitute the specified genes with the phenotype name in the nodes dataframe
            self.nodes['Uniprot'] = self.nodes['Uniprot'].apply(
                lambda x: phenotype_modified if x in unique_uniprot else x)
            self.nodes['Genesymbol'] = self.nodes['Genesymbol'].apply(
                lambda x: phenotype_modified if x in unique_genesymbol else x)

            # Substitute the specified genes with the phenotype name in the edges dataframe
            for column in ['source', 'target']:
                self.edges[column] = self.edges[column].apply(
                    lambda x: phenotype_modified if x in unique_uniprot else x)

            # Group by source and target, and aggregate with the custom function for each column
            self.edges = self.edges.groupby(['source', 'target']).agg({
                'Type': join_unique,  # Aggregate types with the custom function
                'Effect': join_unique,  # Aggregate effects with the custom function
                'References': join_unique  # Aggregate references with the custom function
            }).reset_index()

            # Identify common genes between Uniprot genes and the network's nodes
            common_genes = set(uniprot_genes).intersection(
                set(uniprot_gene_list if uniprot_gene_list else self.nodes["Uniprot"]))
            # For each common gene, add a new edge connecting the gene to the phenotype
            for gene in common_genes:
                new_edge = {"source": gene, "target": phenotype_modified, "Effect": "stimulation",
                            "References": "Gene Ontology"}
                self.edges = self.edges.append(new_edge, ignore_index=True)
        return

    def connect_network_radially(self, max_len: int = 1, direction: Literal['OUT', 'IN', None] = None,
                                 loops: bool = False, consensus: bool = False, only_signed: bool = True) -> None:

        initial_nodes = self.initial_nodes
        initial_nodes_set = set([mapping_node_identifier(i)[2] for i in initial_nodes])

        i = 0
        source_nodes = initial_nodes_set
        target_nodes = initial_nodes_set

        while i < max_len:
            new_nodes = []
            if direction == 'OUT' or direction is None:
                for source in source_nodes:
                    target_neighs = self.connect.find_target_neighbours(source)
                    if source in target_neighs and not loops:
                        target_neighs.remove(source)
                    target_paths = [(source, node) for node in target_neighs]
                    if only_signed:
                        target_paths = self.filter_unsigned_paths(target_paths, consensus)
                    self.add_paths_to_edge_list(target_paths)
                    target_neighs_filtered = [path[1] for path in target_paths]
                    target_neighs_filtered = [node for node in target_neighs_filtered if node not in initial_nodes_set]
                    new_nodes.extend(target_neighs_filtered)
                source_nodes = new_nodes

            new_nodes = []
            if direction == 'IN' or direction is None:
                for target in target_nodes:
                    source_neighs = self.connect.find_source_neighbours(target)
                    if target in source_neighs and not loops:
                        source_neighs.remove(target)
                    source_paths = [(node, target) for node in source_neighs]
                    if only_signed:
                        source_paths = self.filter_unsigned_paths(source_paths, consensus)
                    self.add_paths_to_edge_list(source_paths)
                    source_neighs_filtered = [path[0] for path in source_paths]
                    source_neighs_filtered = [node for node in source_neighs_filtered if node not in initial_nodes_set]
                    new_nodes.extend(source_neighs_filtered)
                target_nodes = new_nodes
            i += 1

        # remove all nodes that have no source or that have no target and are not in the initial nodes
        target_nodes = set(self.edges["target"].unique())
        source_nodes = set(self.edges["source"].unique())

        disconnected_nodes = self.nodes[
            ~self.nodes["Uniprot"].isin(initial_nodes_set) & (
                ~self.nodes["Uniprot"].isin(target_nodes) | ~self.nodes["Uniprot"].isin(source_nodes))]

        # Loop until there are no more disconnected nodes
        while not disconnected_nodes.empty:
            # Collect nodes to remove
            nodes_to_remove = disconnected_nodes["Uniprot"].tolist()

            for node in nodes_to_remove:
                self.remove_node(node)

            # Recalculate disconnected nodes
            target_nodes = set(self.edges["target"].unique())
            source_nodes = set(self.edges["source"].unique())

            disconnected_nodes = self.nodes[
                ~self.nodes["Uniprot"].isin(initial_nodes_set) & (
                    ~self.nodes["Uniprot"].isin(target_nodes) | ~self.nodes["Uniprot"].isin(source_nodes))]

        return

    def connect_as_atopo(self,
                         strategy: Literal['radial', 'complete', None] = None,
                         max_len: int = 1,
                         loops: bool = False,
                         outputs=None,
                         only_signed: bool = True,
                         consensus: bool = False
                         ) -> None:
        """
            This method attempts to connect all nodes of a network object in a topological manner. It iteratively connects
            upstream nodes and checks if the network is connected. If not, it increases the search depth and repeats the process.
            It also removes any nodes that do not have a source in the edge dataframe and are not in the output nodes.

            Parameters:
            - outputs: A list of output nodes to be connected.

            Returns:
            None. The function modifies the network object in-place.
        """

        initial_nodes = [mapping_node_identifier(i)[2] for i in self.initial_nodes]
        initial_nodes_set = set(initial_nodes)

        # Chose the strategy to use to connect the network
        if strategy == 'radial':
            self.connect_network_radially(max_len, direction=None,
                                          loops=loops, consensus=consensus, only_signed=only_signed)
        elif strategy == 'complete':
            self.complete_connection(max_len, minimal=True, k_mean='tight', only_signed=only_signed,
                                     consensus=consensus,
                                     connect_node_when_first_introduced=False)
        else:
            pass

        starting_nodes = set(self.nodes["Uniprot"].tolist())

        if outputs is None:
            return

        # Add output nodes to the network
        for node in outputs:
            self.add_node(node)

        # Convert output nodes to Uniprot identifiers
        outputs_uniprot = [mapping_node_identifier(i)[2] for i in outputs]

        # Initialize depth
        depth = 1

        # While the network is not connected, connect to upstream nodes and first increase the rank and then the depth
        while not is_connected(self):
            self.connect_to_upstream_nodes(outputs_uniprot, depth=depth, rank=len(outputs_uniprot),
                                           only_signed=only_signed, consensus=consensus)
            new_nodes = set(self.nodes["Uniprot"].tolist()) - starting_nodes
            new_nodes = new_nodes - set(outputs_uniprot)

            # Remove nodes that do not have a source in the edge dataframe
            for node in new_nodes:
                if node not in self.edges["target"].unique():
                    self.remove_node(node)
                # remove a node if it auto-regulates itself
                if not loops and any((self.edges['source'] == node) & (self.edges['target'] == node)):
                    self.remove_node(node)

            # If depth reaches 4, stop the process
            if depth == 4:
                print("Current depth is 4, stopping the process")
                break

            # Increase depth
            depth += 1

        # Remove duplicate edges
        self.edges.drop_duplicates().reset_index(drop=True)

        # Create a set of unique sources from the edges DataFrame
        target_nodes = set(self.edges["target"].unique())

        # Identify nodes in the network that are not sources in the edges and not in initial nodes
        disconnected_nodes = self.nodes[
            ~self.nodes["Uniprot"].isin(initial_nodes_set) & ~self.nodes["Uniprot"].isin(target_nodes)]

        # Loop until there are no more disconnected nodes
        while not disconnected_nodes.empty:
            # Collect nodes to remove
            nodes_to_remove = disconnected_nodes["Uniprot"].tolist()

            for node in nodes_to_remove:
                self.remove_node(node)

            # Recalculate disconnected nodes
            target_nodes = set(self.edges["target"].unique())
            disconnected_nodes = self.nodes[
                ~self.nodes["Uniprot"].isin(initial_nodes_set) & ~self.nodes["Uniprot"].isin(target_nodes)]

        return
