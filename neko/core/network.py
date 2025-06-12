from __future__ import annotations
from ..inputs import _universe
from .._methods.enrichment_methods import Connections
import copy
from .._annotations.gene_ontology import Ontology
from .tools import *
from .node import Node
from .edge import Edge
from .network_state import NetworkState
from typing import Optional
from typing_extensions import Literal
from itertools import combinations
from .algorithms.graph_traversal import bfs_algorithm as _bfs_algorithm
from .algorithms.graph_traversal import dfs_algorithm as _dfs_algorithm


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
        self.edges = pd.DataFrame(columns=["source", "target", "Type", "Effect", "References", "Provenance"])
        # Internal object-based storage
        self._node_objs = set()  # Set of Node objects
        self._edge_objs = set()  # Set of Edge objects
        self.initial_nodes = initial_nodes
        self._ontology = Ontology()
        self._populate()
        # --- NetworkState history tracking ---
        self._history: list[NetworkState] = []
        self._history_metadata: list[dict] = []

    def _add_node_obj(self, genesymbol, uniprot, node_type="NaN", metadata=None):
        node_id = uniprot if pd.notna(uniprot) else genesymbol
        node = Node(node_id=node_id, node_type=node_type, metadata=metadata or {})
        self._node_objs.add(node)
        return node

    def _add_edge_obj(self, source, target, interaction_type="undefined", effect=None, references=None, metadata=None):
        edge = Edge(source=source, target=target, interaction_type=interaction_type, evidence=references, metadata=metadata or {})
        self._edge_objs.add(edge)
        return edge

    def nodes_as_objects(self):
        return self._node_objs

    def edges_as_objects(self):
        return self._edge_objs

    def sync_nodes_from_df(self):
        self._node_objs = set()
        for _, row in self.nodes.iterrows():
            self._add_node_obj(row["Genesymbol"], row["Uniprot"], row.get("Type", "NaN"))

    def sync_edges_from_df(self):
        self._edge_objs = set()
        for _, row in self.edges.iterrows():
            self._add_edge_obj(row["source"], row["target"], row.get("Type", "undefined"), row.get("Effect"), row.get("References"))

    def sync_nodes_to_df(self):
        self.nodes = pd.DataFrame([
            {"Genesymbol": n.metadata.get("Genesymbol", n.id), "Uniprot": n.id, "Type": n.type}
            for n in self._node_objs
        ])

    def sync_edges_to_df(self):
        self.edges = pd.DataFrame([
            {"source": e.source, "target": e.target, "Type": e.interaction_type, "Effect": e.metadata.get("Effect"), "References": e.evidence}
            for e in self._edge_objs
        ])

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
                self._add_node_obj(node, node, "NaN")
                return True
            new_entry = {"Genesymbol": genesymbol, "Uniprot": uniprot, "Type": "NaN"}
            self.nodes.loc[len(self.nodes)] = new_entry
            self.nodes = self.nodes.drop_duplicates().reset_index(drop=True)
            self._add_node_obj(genesymbol, uniprot, "NaN")
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
        self._add_node_obj(new_entry["Genesymbol"], new_entry["Uniprot"], new_entry["Type"])
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

    def add_edge(self, edge: pd.DataFrame, provenance: dict = None) -> None:
        """
        This method adds an interaction to the list of interactions while converting it to the NeKo-network format.
        It checks if the edge represents inhibition or stimulation and sets the effect accordingly. It also checks if the
        nodes involved in the interaction are already present in the network, if not, it adds them.
        Optionally, provenance information can be attached to the edge.

        Args:
            - edge: A pandas DataFrame representing the interaction. The DataFrame should contain columns for
            'source', 'target', 'type', and 'references'. The 'source' and 'target' columns represent the nodes involved
            in the interaction. The 'type' column represents the type of interaction. The 'references' column contains
            the references for the interaction.
            - provenance: Optional dict with provenance info (strategy, parameters, timestamp, etc.)

        Returns:
            - None
        """

        # Check if the edge represents inhibition or stimulation and set the effect accordingly
        effect = check_sign(edge)
        references = edge["references"].values[0] if "references" in edge.columns else None
        edge_type = edge["type"].values[0] if "type" in edge.columns else None
        provenance_str = str(provenance) if provenance else None
        df_edge = pd.DataFrame({
            "source": edge["source"],
            "target": edge["target"],
            "Type": edge_type,
            "Effect": effect,
            "References": references,
            "Provenance": provenance_str
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
            # Merge provenance info if needed
            if provenance_str:
                self.edges.loc[existing_edge.index, "Provenance"] = self.edges.loc[existing_edge.index, "Provenance"].astype(str) + "; " + provenance_str
        else:
            # Concatenate the new edge DataFrame with the existing edges in the graph
            self.edges = pd.concat([self.edges, df_edge])

        self.edges = self.edges.drop_duplicates().reset_index(drop=True)
        return

    def get_edge_provenance(self, source: str, target: str) -> str:
        """
        Return the provenance information for a given edge.
        """
        edge_row = self.edges[(self.edges["source"] == source) & (self.edges["target"] == target)]
        if not edge_row.empty:
            return edge_row.iloc[0]["Provenance"]
        return None

    def filter_edges_by_provenance(self, keyword: str) -> pd.DataFrame:
        """
        Return all edges whose provenance contains the given keyword.
        """
        return self.edges[self.edges["Provenance"].str.contains(keyword, na=False)]

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
                "phosphorylate": "undefined",
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

    def connect_nodes(self, only_signed: bool = False, consensus_only: bool = False) -> None:
        """
        Delegates to strategies.connect_nodes.
        """
        from .strategies import connect_nodes
        return connect_nodes(self, only_signed=only_signed, consensus_only=consensus_only)

    def connect_subgroup(self, group, maxlen: int = 1, only_signed: bool = False, consensus: bool = False) -> None:
        """
        Delegates to strategies.connect_subgroup.
        """
        from .strategies import connect_subgroup
        return connect_subgroup(self, group, maxlen=maxlen, only_signed=only_signed, consensus=consensus)

    def connect_component(self, comp_A, comp_B, maxlen: int = 2, mode: Literal['OUT', 'IN', 'ALL'] = 'OUT', only_signed: bool = False, consensus: bool = False) -> None:
        """
        Delegates to strategies.connect_component.
        """
        from .strategies import connect_component
        return connect_component(self, comp_A, comp_B, maxlen=maxlen, mode=mode, only_signed=only_signed, consensus=consensus)

    def connect_to_upstream_nodes(self, nodes_to_connect=None, depth: int = 1, rank: int = 1, only_signed: bool = True, consensus: bool = False) -> None:
        """
        Delegates to strategies.connect_to_upstream_nodes.
        """
        from .strategies import connect_to_upstream_nodes
        return connect_to_upstream_nodes(self, nodes_to_connect=nodes_to_connect, depth=depth, rank=rank, only_signed=only_signed, consensus=consensus)

    def connect_genes_to_phenotype(self, phenotype: str = None, id_accession: str = None, sub_genes: list = None, maxlen: int = 2, only_signed: bool = False, compress: bool = False) -> None:
        """
        Delegates to strategies.connect_genes_to_phenotype.
        """
        from .strategies import connect_genes_to_phenotype
        return connect_genes_to_phenotype(self, phenotype=phenotype, id_accession=id_accession, sub_genes=sub_genes, maxlen=maxlen, only_signed=only_signed, compress=compress)

    def connect_network_radially(self, max_len: int = 1, direction: Literal['OUT', 'IN', None] = None, loops: bool = False, consensus: bool = False, only_signed: bool = True) -> None:
        """
        Delegates to strategies.connect_network_radially.
        """
        from .strategies import connect_network_radially
        return connect_network_radially(self, max_len=max_len, direction=direction, loops=loops, consensus=consensus, only_signed=only_signed)

    def connect_as_atopo(self, strategy: Literal['radial', 'complete', None] = None, max_len: int = 1, loops: bool = False, outputs=None, only_signed: bool = True, consensus: bool = False) -> None:
        """
        Delegates to strategies.connect_as_atopo.
        """
        from .strategies import connect_as_atopo
        return connect_as_atopo(self, strategy=strategy, max_len=max_len, loops=loops, outputs=outputs, only_signed=only_signed, consensus=consensus)

    def complete_connection(self,
                            maxlen: Optional[int] = 2,
                            algorithm: Literal['bfs', 'dfs'] = 'dfs',
                            minimal: bool = True,
                            only_signed: bool = False,
                            consensus: bool = False,
                            connect_with_bias: bool = False,
                            ) -> None:
        """
        Delegates to strategies.complete_connection.
        """
        from .strategies import complete_connection
        return complete_connection(self, maxlen=maxlen, algorithm=algorithm, minimal=minimal, only_signed=only_signed, consensus=consensus, connect_with_bias=connect_with_bias)

    def remove_undefined_interactions(self):
        """
        This function removes all undefined interactions from the network.

        Args:
            - None

        Returns:
            - None
        """
        self.edges = self.edges[self.edges['Effect'] != 'undefined']

    def remove_bimodal_interactions(self):
        """
        Removes all bimodal interactions from the network.
        """
        self.edges = self.edges[self.edges['Effect'] != 'bimodal']

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

    def save_state(self, metadata: Optional[dict] = None) -> None:
        """
        Save the current state of the network (nodes, edges, and optional metadata) to the history.
        """
        state = NetworkState(self.nodes, self.edges, metadata)
        self._history.append(state)
        self._history_metadata.append(state.metadata)

    def restore_state(self, index: int = -1) -> None:
        """
        Restore the network to a previous state by index (default: last state).
        """
        if not self._history:
            raise ValueError("No saved states to restore.")
        state = self._history[index]
        self.nodes = state.nodes.copy(deep=True)
        self.edges = state.edges.copy(deep=True)
        # Optionally, restore other attributes if needed

    def undo(self) -> None:
        """
        Undo the last change by restoring the previous state in history.
        """
        if len(self._history) < 2:
            return
        # Remove the current state
        self._history.pop()
        self._history_metadata.pop()
        # Restore the previous state
        state = self._history[-1]
        self.nodes = state.nodes.copy(deep=True)
        self.edges = state.edges.copy(deep=True)

    def redo(self) -> None:
        """
        Redo is not supported unless branching is implemented. Placeholder for future extension.
        """
        pass

    def compare_states(self, idx1: int, idx2: int) -> dict:
        """
        Compare two saved network states by index. Returns a dict with added/removed nodes/edges and provenance diffs.
        Node comparison is now based only on unique Uniprot identifiers to avoid dtype issues and ensure robust diffing.
        """
        if idx1 >= len(self._history) or idx2 >= len(self._history):
            raise IndexError("State index out of range.")
        state1 = self._history[idx1]
        state2 = self._history[idx2]
        nodes1 = set(str(u) for u in state1.nodes["Uniprot"].drop_duplicates().tolist())
        nodes2 = set(str(u) for u in state2.nodes["Uniprot"].drop_duplicates().tolist())
        edges1 = set(tuple(row) for row in state1.edges.values.tolist())
        edges2 = set(tuple(row) for row in state2.edges.values.tolist())
        diff = {
            "added_nodes": list(nodes2 - nodes1),
            "removed_nodes": list(nodes1 - nodes2),
            "added_edges": list(edges2 - edges1),
            "removed_edges": list(edges1 - edges2),
            "added_provenance": [e for e in (edges2 - edges1) if len(e) > 5 and e[5]],
            "removed_provenance": [e for e in (edges1 - edges2) if len(e) > 5 and e[5]],
        }
        return diff

    def list_states(self) -> list:
        """
        List all saved states (returns their metadata).
        """
        return self._history_metadata.copy()

    def describe_states(self) -> None:
        """
        Print a summary of all saved states.
        """
        for i, meta in enumerate(self._history_metadata):
            print(f"State {i}: {meta}")

    def check_sign(self, edge, consensus_only=False):
        return check_sign(edge, consensus_only)

    def check_gene_list_format(self, genes):
        return check_gene_list_format(genes)

    def mapping_node_identifier(self, node):
        return mapping_node_identifier(node)

    def dfs_algorithm(self,
                      node1: str,
                      node2: str,
                      maxlen: int,
                      only_signed: bool,
                      consensus: bool,
                      connect_with_bias: bool,
                      ) -> None:
        """
        Delegates to graph_traversal.dfs_algorithm.
        """
        _dfs_algorithm(
            find_paths_func=self._connect.find_paths,
            node1=node1,
            node2=node2,
            maxlen=maxlen,
            only_signed=only_signed,
            consensus=consensus,
            connect_with_bias=connect_with_bias,
            add_paths_func=self._add_paths_to_edge_list,
            connect_nodes_func=self.connect_nodes,
            edges_df=self.edges
        )

    def bfs_algorithm(self,
                      node1: str,
                      node2: str,
                      maxlen: int,
                      only_signed: bool,
                      consensus: bool,
                      connect_with_bias: bool,
                      ) -> None:
        """
        Delegates to graph_traversal.bfs_algorithm.
        """
        _bfs_algorithm(
            bfs_func=self._connect.bfs,
            node1=node1,
            node2=node2,
            maxlen=maxlen,
            only_signed=only_signed,
            consensus=consensus,
            connect_with_bias=connect_with_bias,
            add_paths_func=self._add_paths_to_edge_list,
            connect_nodes_func=self.connect_nodes,
            edges_df=self.edges
        )
