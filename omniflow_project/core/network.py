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
import copy
from .._annotations.gene_ontology import Ontology

import pandas as pd


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
        if interaction.get("consensus_stimulation", False):
            return "stimulation"
        elif interaction.get("consensus_inhibition", False):
            return "inhibition"
        else:
            return "undefined"
    else:
        if interaction.get("is_stimulation", False):
            return "stimulation"
        elif interaction.get("is_inhibition", False):
            return "inhibition"
        # Check for "form_complex" column existence
        elif interaction.get("form_complex", False):
            return "form complex"
        else:
            return "undefined"


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
        add_paths_to_edge_list(paths, mode): Adds paths to the edge list of the network.
        connect_nodes(only_signed: bool = False, consensus_only: bool = False): Connects all the nodes in the network.
        is_connected(): Checks if all the nodes in the network are connected.
        filter_unsigned_paths(paths: list[tuple], consensus: bool): Filters out unsigned paths from the provided list of paths.
        connect_subgroup(group: (str | pd.DataFrame | list[str]), maxlen: int = 1, mode: Literal['OUT', 'IN', 'ALL'] = 'ALL', only_signed: bool = False, consensus: bool = False): Connects all the nodes in a particular subgroup.
        complete_connection(maxlen: int = 2, mode: Literal['OUT', 'IN', 'ALL'] = 'ALL', minimal: bool = True, k_mean: Literal['tight', 'extensive'] = 'tight', only_signed: bool = False, consensus: bool = False, connect_node_when_first_introduced: bool = True): Connects all nodes of a network object.
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
            self.resources = resources
        else:
            res = Resources()
            res.load_omnipath_interactions()
            self.resources = res.interactions
        if initial_nodes:
            for node in initial_nodes:
                self.add_node(node)
            self.drop_missing_nodes()
        elif sif_file:
            self.load_network_from_sif(sif_file)

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

    def add_node(self, node: str):
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
        complex_string, genesymbol, uniprot = mapping_node_identifier(node)

        if complex_string:
            new_entry = {"Genesymbol": complex_string, "Uniprot": node, "Type": "NaN"}
        else:
            new_entry = {"Genesymbol": genesymbol, "Uniprot": uniprot, "Type": "NaN"}

        self.nodes.loc[len(self.nodes)] = new_entry
        self.nodes = self.nodes.drop_duplicates()
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

        # add the new nodes to the nodes dataframe
        if not edge["source"].values[0] in self.nodes["Uniprot"].unique():
            self.add_node(edge["source"].values[0])
        if not edge["target"].values[0] in self.nodes["Uniprot"].unique():
            self.add_node(edge["target"].values[0])

        # Concatenate the new edge DataFrame with the existing edges in the graph
        self.edges = pd.concat([self.edges, df_edge])
        self.edges = self.edges.drop_duplicates()
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
                "phosphorilate": "stimulation",
                "stimulation": "stimulation",
                "-1": "inhibition",
                "inhibit": "inhibition",
                "block": "inhibition",
                "inhibition": "inhibition",
                "form complex": "form complex",
                "form_complex": "form complex",
                "form-complex": "form complex",
                "complex formation": "form complex"
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
                    "source": interaction[0],
                    "target": interaction[2],
                    "Type": interaction[3] if len(interaction) > 3 else None,
                    "Effect": effect
                })
                node_set.update([interaction[0], interaction[2]])

        # Create or update the edges DataFrame
        df_edge = pd.DataFrame(interactions)
        self.edges = pd.concat([self.edges, df_edge], ignore_index=True)

        # Update the nodes list
        self.initial_nodes = list(node_set)
        for node in self.initial_nodes:
            self.add_node(node)

        return

    def add_paths_to_edge_list(self, paths, mode):
        """
        Adds paths to the edge list of the network. The paths are added based on the provided mode.

        Parameters:
        - paths: A list of paths, where each path is a sequence of nodes.
        - mode: A string indicating the mode of connection. It can be 'ALL', 'IN', or 'OUT'.

        Returns:
        None. The function modifies the network object in-place.
        """
        database = self.resources

        def add_edge_if_not_empty_and_mode(path, i, mode):
            """
            Helper function to add an edge to the network if the interaction is not empty and matches the provided mode.

            Parameters:
            - path: A sequence of nodes representing a path in the network.
            - i: The index of the current node in the path.
            - mode: A string indicating the mode of connection. It can be 'ALL', 'IN', or 'OUT'.

            Returns:
            None. The function modifies the network object in-place.
            """
            interaction_in = database.loc[(database["source"] == path[i]) &
                                          (database["target"] == path[i + 1])]
            interaction_out = database.loc[(database["target"] == path[i]) &
                                           (database["source"] == path[i + 1])]

            if not interaction_in.empty and (mode in ["ALL", "IN"]):
                self.add_edge(interaction_in)
            elif not interaction_out.empty and (mode in ["ALL", "OUT"]):
                self.add_edge(interaction_out)

        for path in paths:  # Iterate through the list of paths
            if isinstance(path, (str, tuple)):  # Handle single string or tuple
                path = [path]

            for i in range(len(path) - 1):
                add_edge_if_not_empty_and_mode(path, i, mode)

        self.edges = self.edges.drop_duplicates()
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
        self.edges = self.edges.drop_duplicates()

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
            is_full_signed = True
            for i in range(0, len(path)):
                if i == len(path) - 1:
                    break
                interaction_in = interactions.loc[(interactions["source"] == path[i]) &
                                                  (interactions["target"] == path[i + 1])]
                interaction_out = interactions.loc[(interactions["target"] == path[i]) &
                                                   (interactions["source"] == path[i + 1])]
                if not interaction_in.empty and check_sign(interaction_in, consensus) == "undefined":
                    is_full_signed = False
                if not interaction_out.empty and check_sign(interaction_out, consensus) == "undefined":
                    is_full_signed = False
            if is_full_signed:
                filtered_paths.append(path)
        return filtered_paths

    def connect_subgroup(self,
                         group: (str | pd.DataFrame | list[str]),
                         maxlen: int = 1,
                         mode: Literal['OUT', 'IN', 'ALL'] = 'ALL',
                         only_signed: bool = False,
                         consensus: bool = False
                         ):
        """
        This function is used to connect all the nodes in a particular subgroup. It iterates over all pairs of nodes in the
        subgroup and finds paths between them in the resources database. If a path is found, it is added to the edge list of
        the network. The function can operate in different modes ('IN', 'OUT', 'ALL') and can filter out unsigned paths.

        Parameters:
        - group: A list of nodes representing the subgroup to connect. Nodes can be represented as strings, pandas DataFrame, or list of strings.
        - maxlen: The maximum length of the paths to be searched for in the resources database. Default is 1.
        - mode: The mode of connection, which can be 'IN', 'OUT', or 'ALL'. Default is 'ALL'.
        - only_signed: A boolean flag indicating whether to only add signed interactions to the network. Default is False.
        - consensus: A boolean flag indicating whether to only add signed interactions with consensus among references to the network. Default is False.

        Returns:
        None. The function modifies the network object in-place.
        """
        connect = Connections(self.resources)  # here we have the database where we look for the interactions
        if len(self.nodes) == 1:
            print("Number of node insufficient to create connection")
        else:
            for node1, node2 in combinations(group, 2):
                flag = False
                i = 0
                print("looking for paths in the database for node ",
                      mapping_node_identifier(node1)[1], " and ", mapping_node_identifier(node2)[1])
                while not flag and i <= maxlen:

                    if mode == "IN":
                        paths_in = connect.find_paths(node1, node2, maxlen=i,
                                                      mode=mode)
                        paths = paths_in
                    elif mode == "OUT":
                        paths_out = connect.find_paths(node1, node2, maxlen=i,
                                                       mode=mode)
                        paths = paths_out
                    elif mode == "ALL":
                        paths_out = connect.find_paths(node1, node2, maxlen=i,
                                                       mode="OUT")
                        paths_in = connect.find_paths(node1, node2, maxlen=i,
                                                      mode="IN")
                        paths = paths_out + paths_in
                    else:
                        print("The only accepted modes are IN, OUT or ALL, please check the syntax")
                        return
                    if only_signed:
                        paths = self.filter_unsigned_paths(paths, consensus)
                    if not paths:  # if there is no path, look for another one until reach maxlen
                        i += 1
                        continue
                    if paths:
                        print("Found a path!")
                        self.add_paths_to_edge_list(paths, mode)
                        flag = True
        return

    def complete_connection(self,
                            maxlen: int = 2,
                            mode: Literal['OUT', 'IN', 'ALL'] = 'ALL',
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
        - mode: The mode of connection, which can be 'IN', 'OUT', or 'ALL'. Default is 'ALL'.
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
        if k_mean == 'tight':
            i_search = 3
        elif k_mean == 'extensive':
            i_search = 4

        # Create a Connections object for the resources
        connect = Connections(self.resources)

        # Copy the nodes
        nodes = self.nodes.copy()

        # Create a Connections object for the edges
        connect_network = Connections(self.edges)

        # Iterate through all combinations of nodes
        for node1, node2 in combinations(nodes["Uniprot"], 2):
            i = 0
            if minimal:
                # Reset the object connect_network, updating the possible list of paths
                connect_network = Connections(self.edges)
            print("looking for paths in the network for node ",
                  mapping_node_identifier(node1)[1], " and ", mapping_node_identifier(node2)[1])
            while i <= maxlen:

                # As first step, make sure that there is at least one path between two nodes in the network
                if mode == "IN":
                    paths_in = connect_network.find_paths(node1, node2, maxlen=i, mode=mode)
                    paths = paths_in
                elif mode == "OUT":
                    paths_out = connect_network.find_paths(node1, node2, maxlen=i, mode=mode)
                    paths = paths_out
                elif mode == "ALL":
                    paths_out = connect_network.find_paths(node1, node2, maxlen=i, mode="OUT")
                    paths_in = connect_network.find_paths(node1, node2, maxlen=i, mode="IN")
                    paths = paths_out + paths_in
                else:
                    print("The only accepted modes are IN, OUT or ALL, please check the syntax")
                    return

                if paths:
                    print("Found a path!")
                if not paths and i < maxlen:  # if there is no path, look for another one until reach maxlen
                    i += 1
                    continue
                elif not paths and i == maxlen:  # if there is no path of maxlen, look for a new path in the database
                    # until we find a new one
                    flag = False
                    i = 0
                    print("Looking for paths for node ",
                          mapping_node_identifier(node1)[1], " and ", mapping_node_identifier(node2)[1])
                    while not flag and i <= i_search:

                        if mode == "IN":
                            paths_in = connect.find_paths(node1, node2, maxlen=i, mode=mode)
                            paths = paths_in
                        elif mode == "OUT":
                            paths_out = connect.find_paths(node1, node2, maxlen=i, mode=mode)
                            paths = paths_out
                        elif mode == "ALL":
                            paths_out = connect.find_paths(node1, node2, maxlen=i, mode="OUT")
                            paths_in = connect.find_paths(node1, node2, maxlen=i, mode="IN")
                            paths = paths_out + paths_in
                        else:
                            print("The only accepted modes are IN, OUT or ALL, please check the syntax")
                            return
                        if only_signed:
                            paths = self.filter_unsigned_paths(paths, consensus)
                        if not paths:
                            i += 1
                        else:
                            self.add_paths_to_edge_list(paths, mode)
                            if connect_node_when_first_introduced:
                                self.connect_nodes(only_signed, consensus)
                                self.edges = self.edges.drop_duplicates()  # Just in case there are some duplicates
                            flag = True
                    break
                elif paths:  # if there is a path, there is no need to connect the node, so we iterate through another
                    # pair of nodes
                    break
        if not connect_node_when_first_introduced:
            self.connect_nodes(only_signed, consensus)
            self.edges = self.edges.drop_duplicates()  # Just in case there are some duplicates
        return

    def connect_component(self,
                          comp_A: (
                              str | pd.DataFrame | list[str]
                          ),
                          comp_B: (
                              str | pd.DataFrame | list[str]
                          ),
                          maxlen: int = 2,
                          mode: str = 'OUT',
                          only_signed: bool = False,
                          consensus: bool = False):
        """
        This function tries to connect subcomponents of a network object using one of the methods presented in the Connection
        object. This should be a core characteristic of this package and the user should have the possibility to choose
        different methods to enrich its Network object.

        Parameters:
        - comp_A: The first component to connect.
        - comp_B: The second component to connect.
        - maxlen: The maximum length of the paths to be searched for.
        - mode: The mode of connection, which can be 'IN', 'OUT', or 'ALL'.
        - only_signed: A boolean flag to indicate whether to filter unsigned paths.
        - consensus: A boolean flag to indicate whether to check for consensus among references.

        Returns:
        None. The function modifies the network object in-place.
        """
        # Create a Connections object for the resources
        connect = Connections(self.resources)

        # Find paths based on the mode
        if mode in ['IN', 'OUT']:
            paths = connect.find_paths(comp_A, comp_B, maxlen=maxlen, mode=mode)
        elif mode == 'ALL':
            paths = connect.find_paths(comp_A, comp_B, maxlen=maxlen, mode='OUT') + \
                    connect.find_paths(comp_A, comp_B, maxlen=maxlen, mode='IN')
        else:
            print("The only accepted modes are IN, OUT or ALL, please check the syntax")
            return

        # Filter unsigned paths if only_signed is True
        if only_signed:
            paths = self.filter_unsigned_paths(paths, consensus)

        # Add paths to edge list
        self.add_paths_to_edge_list(paths, mode)

        # Calculate the difference between all nodes and comp_A and comp_B
        set_c = list(set(self.nodes['Uniprot'].values) - set(comp_A) - set(comp_B))

        # Connect the subgroup
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
            # Initialize connections object
            connect = Connections(self.resources)
            if nodes_to_connect is None:
                nodes_to_connect = self.nodes["Uniprot"].tolist()

            cascades = connect.find_upstream_cascades(nodes_to_connect, depth, rank)
            if only_signed:
                cascades = self.filter_unsigned_paths(cascades, consensus)
            self.add_cascade_to_edge_list(cascades)
            self.edges.drop_duplicates()
        except Exception as e:
            print(f"An error occurred while connecting to upstream nodes: {e}")
        return

    def convert_edgelist_into_genesymbol(self):
        """
        This function converts the edge dataframe from uniprot to genesymbol.
        """

        def convert_identifier(x):
            identifiers = mapping_node_identifier(x)
            return identifiers[0] or identifiers[1]

        self.edges["source"] = self.edges["source"].apply(convert_identifier)
        self.edges["target"] = self.edges["target"].apply(convert_identifier)

        return

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
        phenotype_genes = self.ontology.get_markers(phenotype=phenotype, id_accession=id_accession)
        if not phenotype_genes:
            print("Something went wrong while getting the markers for: ", phenotype, " and ", id_accession)
            print("Check URL and try again")
            return
        uniprot_genes = [mapping_node_identifier(i)[2] for i in phenotype_genes]

        print("Starting connecting network's nodes to: ", phenotype_genes)
        unique_uniprot = set(uniprot_genes) - set(sub_genes if sub_genes else self.nodes["Uniprot"])
        self.connect_component(sub_genes if sub_genes else self.nodes["Uniprot"].tolist(), list(unique_uniprot),
                               mode="OUT",
                               maxlen=maxlen, only_signed=only_signed)

        if compress:
            phenotype = phenotype or self.ontology.accession_to_phenotype_dict[id_accession]
            phenotype_modified = phenotype.replace(" ", "_")

            # Substitute the specified genes with the phenotype name in the nodes dataframe
            # Substitute the specified genes with the phenotype name in the nodes dataframe
            self.nodes['Uniprot'] = self.nodes['Uniprot'].apply(
                lambda x: phenotype_modified if x in uniprot_genes else x)
            self.nodes['Genesymbol'] = self.nodes['Genesymbol'].apply(
                lambda x: phenotype_modified if x in phenotype_genes else x)

            # Substitute the specified genes with the phenotype name in the edges dataframe
            for column in ['source', 'target']:
                self.edges[column] = self.edges[column].apply(lambda x: phenotype_modified if x in uniprot_genes else x)

            # Group by source and target, and aggregate with the custom function for each column
            self.edges = self.edges.groupby(['source', 'target']).agg({
                'Type': join_unique,  # Aggregate types with the custom function
                'Effect': join_unique,  # Aggregate effects with the custom function
                'References': join_unique  # Aggregate references with the custom function
            }).reset_index()
        return
