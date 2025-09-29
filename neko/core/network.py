from __future__ import annotations
from ..inputs import _universe
from .._methods.enrichment_methods import Connections
import copy
from contextlib import contextmanager
from functools import wraps
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
import pandas as pd
from pandas.util import hash_pandas_object

_METADATA_PREVIEW = 120


def _truncate_repr(value, limit=_METADATA_PREVIEW):
    text = repr(value)
    if len(text) > limit:
        text = text[: limit - 1] + "…"
    return text


def _frame_fingerprint(df: pd.DataFrame) -> tuple:
    if df.empty:
        return (tuple(df.columns), 0, 0)
    hashed = hash_pandas_object(df, index=True)
    return (tuple(df.columns), len(df), int(hashed.sum()))


def _record_state_operation(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        if not getattr(self, "_history_enabled", True):
            return method(self, *args, **kwargs)

        if getattr(self, "_is_initializing", False):
            return method(self, *args, **kwargs)

        depth = getattr(self, "_auto_state_depth", 0)
        self._auto_state_depth = depth + 1

        before_nodes = _frame_fingerprint(self.nodes)
        before_edges = _frame_fingerprint(self.edges)

        try:
            result = method(self, *args, **kwargs)
        finally:
            self._auto_state_depth = depth

        if getattr(self, "_is_initializing", False):
            return result

        if depth > 0:
            return result

        nodes_changed = before_nodes != _frame_fingerprint(self.nodes)
        edges_changed = before_edges != _frame_fingerprint(self.edges)

        if nodes_changed or edges_changed:
            metadata = self._prepare_history_metadata(method.__name__, args, kwargs)
            self.save_state(metadata=metadata)

        return result

    return wrapper


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
        # Internal object-based storage
        self._node_objs = set()  # Set of Node objects
        self._edge_objs = set()  # Set of Edge objects
        self.initial_nodes = initial_nodes
        self._ontology = Ontology()
        # --- NetworkState history tracking ---
        self._states: dict[int, NetworkState] = {}
        self._state_metadata: dict[int, dict] = {}
        self._state_log: list[int] = []
        self._state_counter: int = 0
        self._current_state_id: Optional[int] = None
        self._root_state_id: Optional[int] = None
        self._auto_state_depth: int = 0
        self._history_enabled: bool = True
        self._max_history: Optional[int] = None
        self._is_initializing = True
        self._populate()
        self._is_initializing = False

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
        if self._state_counter == 0:
            self.save_state(metadata={"label": "initial"})

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

    @_record_state_operation
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

    @_record_state_operation
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

    @_record_state_operation
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
        references = edge["references"].values[0] if "references" in edge.columns else None
        edge_type = edge["type"].values[0] if "type" in edge.columns else None
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

    @_record_state_operation
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

    @_record_state_operation
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

    @_record_state_operation
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

    @_record_state_operation
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

    @_record_state_operation
    def connect_nodes(self, only_signed: bool = False, consensus_only: bool = False) -> None:
        """
        Delegates to strategies.connect_nodes.
        """
        from .strategies import connect_nodes
        return connect_nodes(self, only_signed=only_signed, consensus_only=consensus_only)

    @_record_state_operation
    def connect_subgroup(self, group, maxlen: int = 1, only_signed: bool = False, consensus: bool = False) -> None:
        """
        Delegates to strategies.connect_subgroup.
        """
        from .strategies import connect_subgroup
        return connect_subgroup(self, group, maxlen=maxlen, only_signed=only_signed, consensus=consensus)

    @_record_state_operation
    def connect_component(self, comp_A, comp_B, maxlen: int = 2, mode: Literal['OUT', 'IN', 'ALL'] = 'OUT', only_signed: bool = False, consensus: bool = False) -> None:
        """
        Delegates to strategies.connect_component.
        """
        from .strategies import connect_component
        return connect_component(self, comp_A, comp_B, maxlen=maxlen, mode=mode, only_signed=only_signed, consensus=consensus)

    @_record_state_operation
    def connect_to_upstream_nodes(self, nodes_to_connect=None, depth: int = 1, rank: int = 1, only_signed: bool = True, consensus: bool = False) -> None:
        """
        Delegates to strategies.connect_to_upstream_nodes.
        """
        from .strategies import connect_to_upstream_nodes
        return connect_to_upstream_nodes(self, nodes_to_connect=nodes_to_connect, depth=depth, rank=rank, only_signed=only_signed, consensus=consensus)

    @_record_state_operation
    def connect_genes_to_phenotype(self, phenotype: str = None, id_accession: str = None, sub_genes: list = None, maxlen: int = 2, only_signed: bool = False, compress: bool = False) -> None:
        """
        Delegates to strategies.connect_genes_to_phenotype.
        """
        from .strategies import connect_genes_to_phenotype
        return connect_genes_to_phenotype(self, phenotype=phenotype, id_accession=id_accession, sub_genes=sub_genes, maxlen=maxlen, only_signed=only_signed, compress=compress)

    @_record_state_operation
    def connect_network_radially(self, max_len: int = 1, direction: Literal['OUT', 'IN', None] = None, loops: bool = False, consensus: bool = False, only_signed: bool = True) -> None:
        """
        Delegates to strategies.connect_network_radially.
        """
        from .strategies import connect_network_radially
        return connect_network_radially(self, max_len=max_len, direction=direction, loops=loops, consensus=consensus, only_signed=only_signed)

    @_record_state_operation
    def connect_as_atopo(self, strategy: Literal['radial', 'complete', None] = None, max_len: int = 1, loops: bool = False, outputs=None, only_signed: bool = True, consensus: bool = False) -> None:
        """
        Delegates to strategies.connect_as_atopo.
        """
        from .strategies import connect_as_atopo
        return connect_as_atopo(self, strategy=strategy, max_len=max_len, loops=loops, outputs=outputs, only_signed=only_signed, consensus=consensus)

    @_record_state_operation
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

    @_record_state_operation
    def remove_undefined_interactions(self):
        """
        This function removes all undefined interactions from the network.

        Args:
            - None

        Returns:
            - None
        """
        self.edges = self.edges[self.edges['Effect'] != 'undefined']

    @_record_state_operation
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

    def _build_state(self, metadata: Optional[dict], parent_id: Optional[int]) -> NetworkState:
        metadata = metadata or {}
        state = NetworkState(
            nodes=self.nodes,
            edges=self.edges,
            metadata=metadata,
            state_id=self._state_counter,
            parent_ids=[parent_id] if parent_id is not None else [],
        )
        return state

    def save_state(self, metadata: Optional[dict] = None) -> int:
        """Persist the current network snapshot and return its state id."""

        parent_id = self._current_state_id
        state = self._build_state(metadata, parent_id)
        state_id = state.state_id
        self._states[state_id] = state
        self._state_metadata[state_id] = state.metadata
        self._state_log.append(state_id)

        if parent_id is None:
            self._root_state_id = state_id
        else:
            parent_state = self._states[parent_id]
            parent_state.add_child(state_id)

        self._current_state_id = state_id
        self._state_counter += 1
        self._enforce_max_history()
        return state_id

    def _get_state(self, state_id: int) -> NetworkState:
        try:
            return self._states[state_id]
        except KeyError as exc:
            raise ValueError(f"Unknown state id: {state_id}") from exc

    def _resolve_state_id(self, identifier: int) -> int:
        if identifier in self._states:
            return identifier
        if isinstance(identifier, int) and 0 <= identifier < len(self._state_log):
            return self._state_log[identifier]
        raise ValueError(f"Unknown state identifier: {identifier}")

    def set_max_history(self, max_states: Optional[int]) -> None:
        """Set the maximum number of stored history states (None disables pruning)."""

        if max_states is None:
            self._max_history = None
        else:
            max_states = int(max_states)
            if max_states < 2:
                max_states = 2
            self._max_history = max_states
        self._enforce_max_history()

    def _enforce_max_history(self) -> None:
        if self._max_history is None:
            return
        idx = 0
        protected = {self._root_state_id, self._current_state_id}
        while len(self._state_log) > self._max_history and idx < len(self._state_log):
            candidate = self._state_log[idx]
            if candidate in protected:
                idx += 1
                continue
            state = self._states.get(candidate)
            if state is None:
                self._state_log = [sid for sid in self._state_log if sid != candidate]
                idx = 0
                continue
            for child_id in list(state.children_ids):
                child = self._states.get(child_id)
                if child:
                    child.parent_ids = [pid for pid in child.parent_ids if pid != candidate]
                    for parent_id in state.parent_ids:
                        if parent_id not in child.parent_ids:
                            child.parent_ids.append(parent_id)
                        parent_state = self._states.get(parent_id)
                        if parent_state and child_id not in parent_state.children_ids:
                            parent_state.children_ids.append(child_id)
            state.children_ids = []
            for parent_id in state.parent_ids:
                parent = self._states.get(parent_id)
                if parent:
                    parent.children_ids = [c for c in parent.children_ids if c != candidate]
            self._states.pop(candidate, None)
            self._state_metadata.pop(candidate, None)
            self._state_log = [sid for sid in self._state_log if sid != candidate]
            idx = 0

    def set_history_tracking(self, enabled: bool) -> None:
        """Globally enable or disable automatic state capture."""

        self._history_enabled = bool(enabled)

    @contextmanager
    def suspend_history(self):
        """Temporarily suspend automatic state capture within the context."""

        previous = self._history_enabled
        self._history_enabled = False
        try:
            yield
        finally:
            self._history_enabled = previous

    def checkout(self, state_id: int) -> None:
        """Restore the network to a previously saved state."""

        state = self._get_state(state_id)
        self.nodes = state.nodes.copy(deep=True)
        self.edges = state.edges.copy(deep=True)
        self._current_state_id = state_id

    def restore_state(self, state_id: Optional[int] = None) -> None:
        """Backward compatible wrapper around :meth:`checkout`."""

        target = state_id if state_id is not None else self._current_state_id
        if target is None:
            raise ValueError("No saved states to restore.")
        self.checkout(self._resolve_state_id(target))

    def undo(self) -> None:
        """Move to the parent state if available."""

        if self._current_state_id is None:
            return
        parents = self._get_state(self._current_state_id).parent_ids
        if not parents:
            return
        self.checkout(parents[-1])

    def redo(self, state_id: Optional[int] = None) -> None:
        """Move to a child state. If more than one child exists, a specific state id is required."""

        if self._current_state_id is None:
            return
        children = self._get_state(self._current_state_id).children_ids
        if not children:
            return
        target = state_id
        if target is None:
            if len(children) != 1:
                raise ValueError("Multiple branches available; specify a target state id.")
            target = children[0]
        elif target not in children:
            raise ValueError(f"State {target} is not a child of {self._current_state_id}.")
        self.checkout(target)

    def compare_states(self, state_a: int, state_b: int) -> dict:
        """Compare two states identified by their ids."""

        resolved_a = self._resolve_state_id(state_a)
        resolved_b = self._resolve_state_id(state_b)
        state1 = self._get_state(resolved_a)
        state2 = self._get_state(resolved_b)
        nodes1 = self._state_node_labels(state1)
        nodes2 = self._state_node_labels(state2)
        edges1 = self._state_edge_signatures(state1)
        edges2 = self._state_edge_signatures(state2)
        diff = {
            "added_nodes": list(nodes2 - nodes1),
            "removed_nodes": list(nodes1 - nodes2),
            "added_edges": list(edges2 - edges1),
            "removed_edges": list(edges1 - edges2),
        }
        return diff

    def _prepare_history_metadata(self, method_name: str, args, kwargs) -> dict:
        return {
            "method": method_name,
            "args": [self._serialize_history_value(arg) for arg in args],
            "kwargs": {key: self._serialize_history_value(val) for key, val in kwargs.items()},
        }

    def _serialize_history_value(self, value) -> str:
        if isinstance(value, str):
            label = self._node_display_label(value)
            return label if label is not None else value
        if isinstance(value, (list, tuple, set)):
            return "[" + ", ".join(self._serialize_history_value(v) for v in value) + "]"
        if isinstance(value, dict):
            return "{" + ", ".join(f"{k}={self._serialize_history_value(v)}" for k, v in value.items()) + "}"
        return _truncate_repr(value)

    def _state_node_labels(self, state: NetworkState) -> set[str]:
        nodes = state.nodes
        labels: set[str] = set()
        if nodes.empty:
            return labels
        genes = nodes.get("Genesymbol")
        uniprots = nodes.get("Uniprot")
        if genes is not None:
            genes = genes.fillna("")
        if uniprots is not None:
            uniprots = uniprots.fillna("")
        for idx in range(len(nodes)):
            label = ""
            if genes is not None:
                label = str(genes.iloc[idx]) if genes.iloc[idx] else ""
            if (not label) and uniprots is not None:
                label = str(uniprots.iloc[idx]) if uniprots.iloc[idx] else ""
            if not label:
                label = str(nodes.index[idx]) if nodes.index is not None else ""
            if label:
                labels.add(label)
        return labels

    def _node_display_label(self, identifier):
        if pd.isna(identifier):
            return identifier
        identifiers = self.mapping_node_identifier(identifier)
        for candidate in (identifiers[1], identifiers[0], identifiers[2], identifier):
            if candidate:
                return candidate
        return identifier

    def _state_edge_signatures(self, state: NetworkState) -> set[tuple]:
        edges = state.edges
        if edges.empty:
            return set()
        converted = edges.copy()
        if "source" in converted.columns:
            converted["source"] = converted["source"].apply(self._node_display_label)
        if "target" in converted.columns:
            converted["target"] = converted["target"].apply(self._node_display_label)
        return set(tuple(row) for row in converted.to_records(index=False))

    def list_states(self) -> list:
        """Return a creation-ordered list of state metadata."""

        return [
            {"id": state_id, "metadata": self._state_metadata.get(state_id, {})}
            for state_id in self._state_log
        ]

    def describe_history(self) -> None:
        """Pretty-print the branching history tree."""

        if not self._states:
            print("<no states recorded>")
            return

        def _label(state: NetworkState) -> str:
            meta = state.metadata or {}
            label = meta.get("label") or meta.get("description")
            if label:
                return str(label)
            if meta:
                return str(meta)
            return ""

        def _walk(state_id: int, depth: int) -> None:
            state = self._get_state(state_id)
            indent = "  " * depth
            label = _label(state)
            suffix = f" - {label}" if label else ""
            print(f"{indent}State {state_id}{suffix}")
            for child_id in state.children_ids:
                _walk(child_id, depth + 1)

        _walk(self._root_state_id, 0)

    def describe_states(self) -> None:
        """Backward compatible alias for :meth:`describe_history`."""

        self.describe_history()

    @property
    def current_state_id(self) -> Optional[int]:
        return self._current_state_id

    @property
    def root_state_id(self) -> Optional[int]:
        return self._root_state_id

    def history_graph(self):
        """Return a networkx DiGraph representing the state transitions."""

        from .._visual.history import build_history_graph

        return build_history_graph(self)

    def history_digraph(self, include_metadata: bool = True):
        """Return a Graphviz digraph for the history."""

        from .._visual.history import history_digraph

        return history_digraph(self, include_metadata=include_metadata)

    def history_html(self, include_metadata: bool = True, div_class: str = "neko-history-graph") -> str:
        """Return an HTML snippet embedding the history graph."""

        from .._visual.history import history_html

        return history_html(self, include_metadata=include_metadata, div_class=div_class)

    def history_html(self, include_metadata: bool = True, div_class: str = "neko-history-graph") -> str:
        """Return an HTML snippet embedding the history graph."""

        from .._visual.history import history_html

        return history_html(self, include_metadata=include_metadata, div_class=div_class)

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
