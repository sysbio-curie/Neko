from __future__ import annotations
from typing import List, Optional
from itertools import combinations
from ..inputs import _universe
from .._methods.enrichment_methods import Connections
from typing_extensions import Literal
import copy
from .._annotations.gene_ontology import Ontology
from .tools import *



class Network:
    """
    A molecular interaction network.

        The `Network` object is the central organizing component of the `neko`
        module. It is the subject of all operations implemented here, including
        topological algorithms, graph analysis, network visualization and
        integration of database knowledge.

    Args:
        initial_nodes: A list of initial nodes to be added to the network.
        sif_file: A SIF (Simple Interaction Format) file to load the network from.
        resources: A pandas DataFrame containing the resources database.

    Methods:
    """

    def __init__(
            self,
            initial_nodes: list[str] = None,
            sif_file=None,
            resources=None,
        ):

        self._init_args = locals()
        del self._init_args['self']
        self.nodes = pd.DataFrame(columns=["Genesymbol", "Uniprot", "Type"])
        self.edges = pd.DataFrame(columns=["source", "target", "Type", "Effect", "References"])
        self.initial_nodes = initial_nodes
        self._ontology = Ontology()
        self._populate()


    def _populate(self):

        self.resources = (
            _universe.
            network_universe(self._init_args['resources']).
            interactions
        )

        if self.initial_nodes:
            nodes_found = []
            for node in self.initial_nodes:
                if self.add_node(node):
                    nodes_found.append(node)
            self.initial_nodes = nodes_found
            self._drop_missing_nodes()
            self.nodes.reset_index(inplace=True, drop=True)

        elif sif_file := self._init_args['sif_file']:
            self.initial_nodes = []
            self._load_network_from_sif(sif_file)

        self._connect = Connections(self.resources)
        self._algorithms = {
            'dfs': self.dfs_algorithm,
            'bfs': self.bfs_algorithm,
        }

    def copy(self):
        new_instance = copy.deepcopy(self)
        return new_instance

    def check_nodes(self, nodes: list[str]) -> list[str]:
        """
        This function checks if the nodes exist in the resources database and returns the nodes that are present.

        Args:
            - nodes: A list of node identifiers (strings). These are the nodes to be checked.

        Returns:
            - A list[str] of node identifiers that are present in the resources database.
        """
        return [node for node in nodes if
                node in self.resources["source"].unique() or node in self.resources["target"].unique()]

    def check_node(self, node: str) -> bool:
        """
        This function checks if a node exists in the resources' database.

        Args:
            - node: A string representing the node to be checked.

        Returns:
            - A boolean indicating whether the node exists in the resources' database.
        """
        return True if self.check_nodes([node]) else False

    def _drop_missing_nodes(self) -> None:
        """
        This function drops the nodes that are not present in the resources database and print a warning with the
        name of the missing nodes.

        The function works as follows:
        1. It first calls the `check_nodes` function to get a list of nodes that exist in the resources' database.
        2. It then finds the nodes in the network that are not in this list, and removes them from the network.
        3. If there are any missing nodes, it prints a warning with their names.

        This function does not return anything. It modifies the `nodes` attribute of the `Network` object in-place.
        """
        # Get the list of nodes that exist in the resources database
        existing_nodes = list(set(self.check_nodes(self.nodes["Uniprot"].tolist()))
                              |
                              set(self.check_nodes(self.nodes["Genesymbol"].tolist())))

        # Create a mask for nodes that exist in either Uniprot or Genesymbol lists
        existing_mask = self.nodes["Uniprot"].isin(existing_nodes) | self.nodes["Genesymbol"].isin(existing_nodes)
        # Identify removed nodes
        removed_nodes = self.nodes[~existing_mask]

        # Print a warning with the name of the missing nodes
        if not removed_nodes.empty:
            missing_uniprot = removed_nodes["Uniprot"].tolist()
            missing_genesymbol = removed_nodes["Genesymbol"].tolist()
            missing_identifiers = list(set(missing_uniprot + missing_genesymbol))  # Remove duplicates
            print(
                "Warning: The following nodes were not found in the resources database and have been removed from the "
                "network:",
                ", ".join(str(node) for node in missing_identifiers if pd.notna(node))
            )
        # Keep only the nodes that exist in the resources database
        self.nodes = self.nodes[existing_mask]

        return

    def add_node(self, node: str, from_sif: bool = False) -> bool:
        """
        Adds a node to the network. The node is added to the nodes DataFrame of the network. The function checks the
        syntax for the genesymbol to ensure it is correct. If the node is a complex, it is added with the
        'Genesymbol' as the complex string and 'Uniprot' as the node. Otherwise, it is added with the 'Genesymbol' as
        the genesymbol and 'Uniprot' as the uniprot. The 'Type' is set as 'NaN' for all new nodes.

        Args:
            - node: A string representing the node to be added. The node can be represented by either its
                    Genesymbol or Uniprot identifier.

        Returns:
            - None.
        """

        if from_sif:
            # check that the new entry node can be translated using the function mapping node identifier (all the
            # output of the function should be None) if it cannot be translated, print an error message but add the
            # node to the network anyway

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
            return True

        complex_string, genesymbol, uniprot = mapping_node_identifier(node)

        if complex_string:
            new_entry = {"Genesymbol": complex_string, "Uniprot": node, "Type": "NaN"}
        else:
            new_entry = {"Genesymbol": genesymbol, "Uniprot": uniprot, "Type": "NaN"}

        if not self.check_node(uniprot) and not self.check_node(genesymbol):
            print("Error: node %s is not present in the resources database" % node)
            return False

        self.nodes.loc[len(self.nodes)] = new_entry

        self.nodes = self.nodes.drop_duplicates().reset_index(drop=True)
        return True

    def remove_node(self, node: str) -> None:
        """
        Removes a node from the network. The node is removed from both the list of nodes and the list of edges.

        Args:
            - node: A string representing the node to be removed. The node can be represented by either its
                    Genesymbol or Uniprot identifier.

        Returns:
            - None
        """
        # Remove the node from the nodes DataFrame
        self.nodes = self.nodes[(self.nodes.Genesymbol != node) & (self.nodes.Uniprot != node)]

        # Translate the node identifier to Uniprot
        node = mapping_node_identifier(node)[2]

        # Remove any edges associated with the node from the edges DataFrame
        self.edges = self.edges[~self.edges[['source', 'target']].isin([node]).any(axis=1)]

        return

    def add_edge(self, edge: pd.DataFrame) -> None:
        """
        This method adds an interaction to the list of interactions while converting it to the NeKo-network format.
        It checks if the edge represents inhibition or stimulation and sets the effect accordingly. It also checks if the
        nodes involved in the interaction are already present in the network, if not, it adds them.

        Args:
            - edge: A pandas DataFrame representing the interaction. The DataFrame should contain columns for
            'source', 'target', 'type', and 'references'. The 'source' and 'target' columns represent the nodes involved
            in the interaction. The 'type' column represents the type of interaction. The 'references' column contains
            the references for the interaction.

        Returns:
            - None
        """

        # Check if the edge represents inhibition or stimulation and set the effect accordingly
        effect = check_sign(edge)
        # check if the column reference is present in the edge dataframe
        references = edge["references"].values[0] if "references" in edge.columns else None
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
        if not existing_edge.empty and references is not None:
            self.edges.loc[existing_edge.index, "References"] += "; " + str(references)
        else:
            # Concatenate the new edge DataFrame with the existing edges in the graph
            self.edges = pd.concat([self.edges, df_edge])

        self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def remove_edge(self, node1: str, node2: str) -> None:
        """
        This function removes an edge from the network. It takes the source node and target node as input and removes
        the edge from the edges DataFrame.

        Args:
            - node1: A string representing the source node of the edge.
            - node2: A string representing the target node of the edge.

        Returns:
            - None
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

    def remove_disconnected_nodes(self) -> None:
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
                         type: Literal['Genesymbol', 'Uniprot', 'both'] = 'Genesymbol'
                         ) -> None:
        """
        This function modifies the name of a node in the network. It takes the old name of the node and the new name
        as input and modifies the name of the node in the nodes and in the edges DataFrame. If type is set to
        'Genesymbol', it modifies the genesymbol name of the node in the nodes DataFrame. If type is set to
        'Uniprot', it modifies the uniprot name of the node in the edges DataFrame. If type is set to 'both',
        it modifies both the genesymbol and uniprot names of the node in the nodes and edges DataFrame.


        Args:
            - old_name: A string representing the old name of the node. - new_name: A string representing the new
            name of the node. - type: A string indicating the type of name to be modified. It can be 'Genesymbol',
            'Uniprot', or 'both'. Default is 'Genesymbol'.

        Returns:
            -None
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

    def print_my_paths(self,
                       node1: str,
                       node2: str,
                       maxlen: int = 2,
                       genesymbol: bool = True
                       ) -> None:
        """
        This function prints all the paths between two nodes in the network. It uses the `find_paths` method from the
        `Connections` class to find all the paths between the two nodes in the Network object. If no paths are found,
        it prints a warning message. If one of the selected nodes is not present in the network, it prints an error message.

        Args:
            - node1: A string representing the source node.
            - node2: A string representing the target node.
            - maxlen: An integer representing the maximum length of the paths to be searched for. Default is 2.
            - genesymbol: A boolean flag indicating whether to print the paths in genesymbol format. Default is True.

        Returns:
            - None
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

    def remove_path(self, path: list[str]) -> None:
        """
        This function removes a path from the network. It takes a list of nodes representing the path and removes all
        the edges between the nodes in the path.

        Args:
            - path: A list of nodes representing the path to be removed. The nodes can be represented as strings or tuples.

        Returns:
            - None
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

    def _load_network_from_sif(self, sif_file) -> None:
        """
        Load a network object from a SIF (Simple Interaction Format) file.

        Args:
            - sif_file: A string representing the path to the SIF file.

        Returns:
            None
        """
        interactions = []
        node_set = set()

        def determine_effect(interaction_type):
            """
            Determine the effect based on the interaction type.

            Args:
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
                "complex_formation": "form complex",
                "bimodal": "bimodal",
                "both": "bimodal"
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

    def _add_paths_to_edge_list(self, paths) -> None:
        """
        This method adds paths to the edge list of the network. A path is a sequence of nodes where each node is
        connected to the next node in the sequence. The function checks if there is an interaction between each pair
        of nodes in the path in the resources' database. If an interaction exists, it is added to the edge list of
        the network.

        Args:
            - paths: A list of paths, where each path is a sequence of nodes. A node can be a string or a tuple.

        Returns:
            - None
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
                    if not ((self.edges['source'] == interaction['source'].values[0]) &
                            (self.edges['target'] == interaction['target'].values[0])).any():
                        self.add_edge(interaction)

        # Remove duplicate edges from the edge list
        self.edges = self.edges.drop_duplicates().reset_index(drop=True)

        return

    def _add_cascade_to_edge_list(self, cascades) -> None:
        """
        This function adds cascades to the edge list of the network. A cascade is a sequence of nodes where each node is
        connected to the next node in the sequence. The function checks if there is an interaction between each pair of nodes
        in the cascade in the resources' database. If an interaction exists, it is added to the edge list of the network.

        Args:
            - cascades: A list of cascades, where each cascade is a sequence of nodes.

        Returns:
            - None
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
                      consensus_only: bool = False
                      ) -> None:
        """
        Basic node connections. It adds all the interactions found in the omnipath database.
        Once an interaction is found it will be added to the list of edges.
        The only_signed flag makes sure that just signed interaction will be added to the network, while "consensus_only"
        makes sure that just signed interaction with consensus among references will be included.

        Args:
            - only_signed: A boolean flag indicating whether to only add signed interactions to the network.
            - consensus_only: A boolean flag indicating whether to only add signed interactions with consensus among
                            references to the network.

        Returns:
            - None
        """

        if len(self.nodes) == 1:
            print("Number of node insufficient to create connection")
            return

        def add_edge_if_not_empty_and_signed(node1, node2):
            """
            Helper function to add an edge to the network if the interaction is not empty and, if the `only_signed`
            flag is set, the interaction is signed.

            Args:
            - node1: The source node of the interaction.
            - node2: The target node of the interaction.

            Returns:
            None. The function modifies the network object in-place.
            """
            if node2 in self._connect.find_all_neighbours(node1):
                interaction = self.resources.loc[(self.resources["source"] == node1) &
                                                 (self.resources["target"] == node2)]
                if not interaction.empty and (
                    not only_signed or check_sign(interaction, consensus_only) != "undefined"):
                    self.add_edge(interaction)

        for node1, node2 in combinations(self.nodes["Uniprot"], 2):
            add_edge_if_not_empty_and_signed(node1, node2)
            add_edge_if_not_empty_and_signed(node2, node1)
        return

    def is_connected(self) -> bool:
        """
        This function checks if the network is connected. It uses Depth-First Search (DFS) to traverse the network
        and check if all nodes are visited. If all nodes are visited, the network is connected.

        Args:
            - None

        Returns:
            - None

        """

        # Create a set to store visited nodes
        visited = set()

        # Function for Depth-First Search
        def dfs(node):
            visited.add(node)
            for neighbour in self._connect.find_all_neighbours(node):
                if neighbour not in visited:
                    dfs(neighbour)

        # Start DFS from the first node
        dfs(self.nodes.iloc[0]['Uniprot'])

        # Check if all nodes are visited
        return set(self.nodes['Uniprot']) == visited

    def remove_undefined_interactions(self):
        """
        This function removes all undefined interactions from the network.

        Args:
            - None

        Returns:
            - None
        """
        self.edges = self.edges[self.edges['Effect'] != 'undefined']

    def _filter_unsigned_paths(self,
                                paths: list[list[str]],
                                consensus: bool
                                ) -> list[list[str]]:
        """
        This function filters out unsigned paths from the provided list of paths. An unsigned path is a path where at
        least one interaction does not have a defined sign (stimulation or inhibition). The function checks each
        interaction in each path, and if the interaction is unsigned, the path is not included in the output list.

        Args:
            - paths: A list of paths, where each path is a sequence of nodes.
            - consensus: A boolean indicating whether to check for consensus among references when determining the sign of an interaction.

        Returns:
            - A list[tuple] of paths where all interactions in each path are signed.
        """

        interactions = self.resources
        filtered_paths = []
        for path in paths:
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

    def convert_edgelist_into_genesymbol(self) -> pd.DataFrame:
        """
        This function generates a new edges dataframe with the source and target identifiers translated (if possible)
        in Genesymbol format.

        Args:
             - None

        Returns:
            - A pandas DataFrame containing the edges with the source and target identifiers translated into Genesymbol
                format.
        """

        def convert_identifier(x):
            identifiers = mapping_node_identifier(x)
            return identifiers[0] or identifiers[1]

        gs_edges = self.edges.copy()

        gs_edges["source"] = gs_edges["source"].apply(convert_identifier)
        gs_edges["target"] = gs_edges["target"].apply(convert_identifier)

        return gs_edges

    ###################################################################################################################
    ######################################### STRATEGY CONNECTIONS ####################################################
    ###################################################################################################################

    def connect_subgroup(self,
                         group: (str | pd.DataFrame | list[str]),
                         maxlen: int = 1,
                         only_signed: bool = False,
                         consensus: bool = False
                         ) -> None:
        """
        This function is used to connect all the nodes in a particular subgroup. It iterates over all pairs of nodes
        in the subgroup and finds paths between them in the resources' database. If a path is found, it is added to
        the edge list of the network. The function also filters out unsigned paths if the `only_signed` flag is set
        to True.

        Args:
            - group: A list of nodes representing the subgroup to connect. Nodes can be represented as strings,
                        pandas DataFrame, or list of strings.
            - maxlen: The maximum length of the paths to be searched for in the resources' database. Default is 1.
            - only_signed: A boolean flag indicating whether to only add signed interactions to the network. Default is
                            False.
            - consensus: A boolean flag indicating whether to only add signed interactions with consensus among
                            references to the network. Default is False.

        Returns:
            - None
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
                        paths_out = self._connect.find_paths(node1, node2, maxlen=i)
                        if only_signed:
                            paths_out = self._filter_unsigned_paths(paths_out, consensus)
                    if not paths_in:
                        paths_in = self._connect.find_paths(node2, node1, maxlen=i)
                        if only_signed:
                            paths_in = self._filter_unsigned_paths(paths_in, consensus)
                    if not paths_in or not paths_out and i <= maxlen:
                        i += 1
                    if (paths_in or paths_out) and i > maxlen or (paths_in and paths_out):
                        paths = paths_out + paths_in
                        self._add_paths_to_edge_list(paths)
                        break
        return

    def dfs_algorithm(self,
                      node1: str,
                      node2: str,
                      maxlen: int,
                      only_signed: bool,
                      consensus: bool,
                      connect_with_bias: bool,
                      ) -> None:

        """

        This function uses the Depth-First Search (DFS) algorithm to find paths between two nodes in the network. It
        starts from the target node and searches for paths of increasing length until it finds a path of length
        maxlen. If the `only_signed` flag is set to True, it filters out unsigned paths. If the `connect_with_bias`
        flag is set to True, it connects the nodes when first introduced. Args: node1: node2: maxlen: only_signed:
        consensus: connect_with_bias:

        Args:
            - node1: A string representing the source node.
            - node2: A string representing the target node.
            - maxlen: An integer representing the maximum length of the paths to be searched for. Default is 2.
            - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
            - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
            - connect_with_bias: A boolean flag indicating whether to connect the nodes when first introduced.
                                Default is True.

        Returns:
            - None
        """


        paths = self._connect.find_paths(start=node1, end=node2, maxlen=maxlen, minlen=1)
        if only_signed:
            paths = self._filter_unsigned_paths(paths, consensus)
        if paths:
            self._add_paths_to_edge_list(paths)
            if connect_with_bias:
                self.connect_nodes(only_signed, consensus)
                self.edges = self.edges.drop_duplicates().reset_index(drop=True)

    def bfs_algorithm(self,
                      node1: str,
                      node2: str,
                      maxlen: int,
                      only_signed: bool,
                      consensus: bool,
                      connect_with_bias: bool,
                      ) -> None:

        """

        This function uses the Breadth-First Search (BFS) algorithm to find paths between two nodes in the network.
        It starts from the target node and searches for paths of increasing length until it finds a path of length
        maxlen. If the `only_signed` flag is set to True, it filters out unsigned paths. If the `connect_with_bias`
        flag is set to True, it connects the nodes when first introduced. Args: node1: node2: only_signed: consensus:
        connect_with_bias:

        Args:
            - node1: A string representing the source node.
            - node2: A string representing the target node.
            - maxlen: An integer representing the maximum length of the paths to be searched for. Default is 2.
            - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
            - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
            - connect_with_bias: A boolean flag indicating whether to connect the nodes when first introduced.
                                Default is True.

        Returns:
            - None

        """

        paths = self._connect.bfs(start=node1, end=node2, maxlen=maxlen)
        if only_signed:
            paths = self._filter_unsigned_paths(paths, consensus)
        if paths:
            self._add_paths_to_edge_list(paths)
            if connect_with_bias:
                self.connect_nodes(only_signed, consensus)
                self.edges = self.edges.drop_duplicates().reset_index(drop=True)



    def complete_connection(self,
                            maxlen: Optional[int] = 2,
                            algorithm: Literal['bfs', 'dfs'] = 'dfs',
                            minimal: bool = True,
                            only_signed: bool = False,
                            consensus: bool = False,
                            connect_with_bias: bool = False,
                            ) -> None:
        """
        This function attempts to connect all nodes of a network object using one of the methods presented in the
        Connection object. This is a core characteristic of this package and the user should have the possibility to
        choose different methods to enrich its Network object.

        Args:
            - maxlen: The maximum length of the paths to be searched for. Default is 2.
            - algorithm: The search algorithm to be used. It can be 'bfs' (Breadth-First Search) or 'dfs'
                (Depth-First Search).
            - minimal: A boolean flag indicating whether to reset the object connect_network, updating the possible list
                of paths. Default is True.
            - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
            - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
            - connect_with_bias: A boolean flag indicating whether to connect nodes when first
                introduced. Default is True.

        Returns:
            - None
        """

        # Copy the nodes
        nodes = self.nodes.copy()

        # Create a Connections object for the edges
        connect_network = Connections(self.edges)

        # Iterate through all combinations of nodes
        for node1, node2 in combinations(nodes["Uniprot"], 2):
            if not self.check_node(node1) or not self.check_node(node2):
                print(
                    "Error: node %s is not present in the resources database" % node1 if not self.check_node(
                        node1) else "Error: node %s is not present in the resources database" % node2)
                continue
            i = 0
            # Reset the object connect_network, updating the possible list of paths if minimal is True
            if minimal:
                connect_network = Connections(self.edges)

            # As first step, make sure that there is at least one path between two nodes in the network
            paths_in = connect_network.bfs(start=node2, end=node1, maxlen=maxlen)
            paths_out = connect_network.bfs(start=node1, end=node2, maxlen=maxlen)

            if not paths_in:
                self._algorithms[algorithm](node1=node2, node2=node1, maxlen=maxlen, only_signed=only_signed, consensus=consensus, connect_with_bias=connect_with_bias)
            if not paths_out:
                self._algorithms[algorithm](node1=node1, node2=node2, maxlen=maxlen, only_signed=only_signed, consensus=consensus, connect_with_bias=connect_with_bias)

        # If connect_with_bias is False, connect nodes after all paths have been found
        if not connect_with_bias:
            self.connect_nodes(only_signed, consensus)
            self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def connect_component(self,
                          comp_A: (str | list[str]),
                          comp_B: (str | list[str]),
                          maxlen: int = 2,
                          mode: Literal['OUT', 'IN', 'ALL'] = 'OUT',
                          only_signed: bool = False,
                          consensus: bool = False
                          ) -> None:
        """
        This function attempts to connect subcomponents of a network object using one of the methods presented in the
        Connection object. This is a core characteristic of this package and the user should have the possibility to
        choose different methods to enrich its Network object.

        Args:
            - comp_A: A string or list of strings representing the first component to connect.
            - comp_B: A string or list of strings representing the second component to connect.
            - maxlen: The maximum length of the paths to be searched for. Default is 2.
            - mode: The search mode, which can be 'OUT', 'IN', or 'ALL'. Default is 'OUT'.
            - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
            - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.

        Returns:
            - None

        """

        # Determine the search mode and find paths accordingly
        if mode == "IN":
            paths_in = self._connect.find_paths(comp_B, comp_A, maxlen=maxlen)
            paths = paths_in
        elif mode == "OUT":
            paths_out = self._connect.find_paths(comp_A, comp_B, maxlen=maxlen)
            paths = paths_out
        elif mode == "ALL":
            paths_out = self._connect.find_paths(comp_A, comp_B, maxlen=maxlen)
            paths_in = self._connect.find_paths(comp_B, comp_A, maxlen=maxlen)
            paths = paths_out + paths_in
        else:
            print("The only accepted modes are IN, OUT or ALL, please check the syntax")
            return

        # Filter unsigned paths if the only_signed flag is set
        if only_signed:
            paths = self._filter_unsigned_paths(paths, consensus)

        # Add the paths to the edge list
        self._add_paths_to_edge_list(paths)

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
                                  consensus: bool = False
                                  ) -> None:
        """
        This function connects the provided nodes to their upstream nodes in the network.

        Args:
            - nodes_to_connect: A list of nodes to connect. If not provided, all nodes in the network are considered.
            - depth: The depth of the search for upstream nodes.
            - rank: The rank of the search for upstream nodes.
            - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is True.
            - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.

        Returns:
            - None
        """
        try:
            if nodes_to_connect is None:
                nodes_to_connect = self.nodes["Uniprot"].tolist()

            cascades = self._connect.find_upstream_cascades(nodes_to_connect, depth, rank)

            if only_signed:
                cascades = self._filter_unsigned_paths(cascades, consensus)
            self._add_cascade_to_edge_list(cascades)
            self.edges.drop_duplicates().reset_index(drop=True)
        except Exception as e:
            print(f"An error occurred while connecting to upstream nodes: {e}")
        return

    def connect_genes_to_phenotype(self,
                                   phenotype: str = None,
                                   id_accession: str = None,
                                   sub_genes: list[str] = None,
                                   maxlen: int = 2,
                                   only_signed: bool = False,
                                   compress: bool = False
                                   ) -> None:
        """
        This function connects genes to a phenotype based on the provided Args. It retrieves phenotype markers,
        identifies unique Uniprot genes, and connects them to the network. It also has the option to compress the
        network by substituting specified genes with the phenotype name.

        Args:
            - phenotype: The phenotype to connect to. If not provided, it will be retrieved using the id_accession.
            - id_accession: The accession id of the phenotype. If not provided, the phenotype parameter must be given.
            - sub_genes: A list of genes to be considered for connection. If not provided, all nodes in the network are
                considered.
            - maxlen: The maximum length of the paths to be searched for.
            - only_signed: A boolean flag to indicate whether to filter unsigned paths.
            - compress: A boolean flag to indicate whether to substitute the specified genes with the phenotype name.

        Returns:
            - None
        """
        # Initialize lists for Uniprot and genesymbol genes
        uniprot_gene_list = []
        genesymbols_genes = []

        # Retrieve phenotype markers
        phenotype_genes = self._ontology.get_markers(phenotype=phenotype, id_accession=id_accession)
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
            phenotype = phenotype or self._ontology.accession_to_phenotype_dict[id_accession]
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
                'Effect': determine_most_frequent_effect,  # Use the new function to determine the most frequent effect
                'References': join_unique  # Aggregate references with the custom function
            }).reset_index()

            # Identify common genes between Uniprot genes and the network's nodes
            common_genes = set(uniprot_genes).intersection(
                set(uniprot_gene_list if uniprot_gene_list else self.nodes["Uniprot"]))
            # For each common gene, add a new edge connecting the gene to the phenotype
            for gene in common_genes:
                new_edge = pd.DataFrame({"source": [gene], "target": [phenotype_modified], "Effect": ["stimulation"],
                                         "References": ["Gene Ontology"]})
                self.edges = pd.concat([self.edges, new_edge], ignore_index=True)
        return

    def connect_network_radially(self,
                                 max_len: int = 1,
                                 direction: Literal['OUT', 'IN', None] = None,
                                 loops: bool = False,
                                 consensus: bool = False,
                                 only_signed: bool = True
                                 ) -> None:

        """
        This function connects all nodes of a network object in a radial manner. It iteratively connects upstream and
        downstream nodes of the initial nodes. The function also removes any nodes that do not have a source in the edge
        dataframe and are not in the initial nodes.

        Args:
            - max_len: The maximum length of the paths to be searched for. Default is 1.
            - direction: The direction of the search. It can be 'OUT', 'IN', or None. Default is None.
            - loops: A boolean flag indicating whether to allow loops in the network. Default is False.
            - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.
            - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is True.

        Returns:
            - None

        """

        initial_nodes = self.initial_nodes
        initial_nodes_set = set([mapping_node_identifier(i)[2] for i in initial_nodes])

        i = 0
        source_nodes = initial_nodes_set
        target_nodes = initial_nodes_set

        while i < max_len:
            new_nodes = []
            if direction == 'OUT' or direction is None:
                for source in source_nodes:
                    target_neighs = self._connect.find_target_neighbours(source)
                    if source in target_neighs and not loops:
                        target_neighs.remove(source)
                    target_paths = [(source, node) for node in target_neighs]
                    if only_signed:
                        target_paths = self._filter_unsigned_paths(target_paths, consensus)
                    self._add_paths_to_edge_list(target_paths)
                    target_neighs_filtered = [path[1] for path in target_paths]
                    target_neighs_filtered = [node for node in target_neighs_filtered if node not in initial_nodes_set]
                    new_nodes.extend(target_neighs_filtered)
                source_nodes = new_nodes

            new_nodes = []
            if direction == 'IN' or direction is None:
                for target in target_nodes:
                    source_neighs = self._connect.find_source_neighbours(target)
                    if target in source_neighs and not loops:
                        source_neighs.remove(target)
                    source_paths = [(node, target) for node in source_neighs]
                    if only_signed:
                        source_paths = self._filter_unsigned_paths(source_paths, consensus)
                    self._add_paths_to_edge_list(source_paths)
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
        This method attempts to connect all nodes of a network object in a topological manner. It iteratively
        connects upstream nodes and checks if the network is connected. If not, it increases the search depth and
        repeats the process. It also removes any nodes that do not have a source in the edge dataframe and are not in
        the output nodes.

        Args:
            - strategy: The strategy to use to connect the network. It can be 'radial' or 'complete'. Default is None.
            - max_len: The maximum length of the paths to be searched for. Default is 1.
            - loops: A boolean flag indicating whether to allow loops in the network. Default is False.
            - outputs: A list of output nodes to connect to. Default is None.
            - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is True.
            - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.

        Returns:
            - None
        """

        initial_nodes = [mapping_node_identifier(i)[2] for i in self.initial_nodes]
        initial_nodes_set = set(initial_nodes)

        # Chose the strategy to use to connect the network
        if strategy == 'radial':
            self.connect_network_radially(max_len, direction=None,
                                          loops=loops, consensus=consensus, only_signed=only_signed)
        elif strategy == 'complete':
            self.complete_connection(max_len, minimal=True, only_signed=only_signed,
                                     consensus=consensus,
                                     connect_with_bias=False)
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
