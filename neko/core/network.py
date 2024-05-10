from __future__ import annotations
from typing import List, Optional, Tuple
from itertools import combinations
from typing_extensions import Literal
import copy

from tqdm import tqdm
import pandas as pd
from pypath.utils import mapping
import omnipath as op

from .._inputs._db.omnipath import Resources
from .._methods.enrichment_methods import Connections
from .._annotations import go
import _networkbase as _nbase


class Network(_nbase.NetworkBase):
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

    def __init__(
            self,
            initial_nodes: list[str] = None,
            sif_file=None,
            resources=None
    ):
        self.nodes = pd.DataFrame(columns=["Genesymbol", "Uniprot", "Type"])
        self.edges = pd.DataFrame(columns=["source", "target", "Type", "Effect", "References"])
        self.initial_nodes = initial_nodes
        self.ontology = None
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
        return self.resources & set(nodes)


    def __iadd__(self, other):

        if isinstance(other, str):

            self.add_node(other)

        else:

            raise ValueError('We can add only strings (node IDs).')


    def __isub__(self, other):

        if isinstance(other, str):

            self.remove_node(other)

        else:

            raise ValueError('We can remove only strings (node IDs).')


    def __add__(self, other) -> Network:

        self.__iadd__(other)
        return self


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
        complex_string, genesymbol, uniprot = translate_id(node)

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
        node = translate_id(node)[2]

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
            node1 = translate_id(node1)[2]
        if check_gene_list_format([node2]):
            node2 = translate_id(node2)[2]

        # Remove the edge from the edges DataFrame, if the effect or the nodes are not present, print a warning
        if not self.edges[(self.edges["source"] == node1) & (self.edges["target"] == node2)].empty:
            self.edges = self.edges[~((self.edges["source"] == node1) & (self.edges["target"] == node2))]
        else:
            print("Warning: The edge does not exist in the network, check syntax for ",
                  translate_id(node1)[1], " and ", translate_id(node2)[1])
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
            node1 = translate_id(node1)[2]
        if check_gene_list_format([node2]):
            node2 = translate_id(node2)[2]

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

    @classmethod
    def from_sif(cls, path: str):
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
                "->": "stimulation",
                "-|": "inhibition",
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
                    "source": translate_id(interaction[0])[2] if check_gene_list_format(
                        [interaction[0]]) else interaction[0],
                    "target": translate_id(interaction[2])[2] if check_gene_list_format(
                        [interaction[2]]) else interaction[2],
                    "Type": interaction[3] if len(interaction) > 3 else None,
                    "Effect": effect,
                    "References": "SIF file"
                })
                node_set.update([interaction[0], interaction[2]])

        # Create or update the edges DataFrame
        df_edge = pd.DataFrame(interactions)

        new = cls()
        new.edges = pd.concat([self.edges, df_edge], ignore_index=True)
        new.initial_nodes = list(node_set)
        for node in new.initial_nodes:
            new.add_node(node)

        return new


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


    def check_node_existence(self, node: str) -> bool:
        """
        This function checks if a node exists in the resources' database.

        Parameters:
        - node: A string representing the node to be checked.

        Returns:
        - A boolean indicating whether the node exists in the resources' database.
        """
        return node in self.resources["source"].unique() or node in self.resources["target"].unique()
