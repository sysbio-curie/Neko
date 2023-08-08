from pypath.utils import mapping
import omnipath as op
from itertools import combinations
import pandas as pd
from .._inputs.resources import Resources
from .._methods.enrichment_methods import Connections
from typing_extensions import Literal


class Network:
    """
    Main class of the Omniflow package, it stores nodes and edges and offers different methods for enrichment analysis.
    It should take as main argument a list of nodes and as optional a series of filters and others options (inputs/outputs,
    specific tissues, type of interactions, prune inputs/outputs etc...).
    For the enrichment part, the user should select one or more databases to retrieve the desired interactions. Those
    databases come from the omnipath/pypath module or from a local database provided by the user.
    Nodes are stored in a list Node object.
    Edges are stored in a list of Edge objects.
    """

    def __init__(self,
                 initial_nodes: list[str]):
        self.nodes = pd.DataFrame(columns=["Genesymbol", "Uniprot", "Type"])
        self.edges = pd.DataFrame(columns=["source", "target", "Type", "Effect"])
        res = Resources()
        res.load_omnipath_interactions()
        self.resources = res.omnipath_interactions  ### in future release, a string can determine which database to use/load
        if initial_nodes:
            for node in initial_nodes:
                self.add_node(node)

    def add_node(self,
                 node: str):
        """
        Add a node to the list of nodes checking that the syntax for the genesymbol is actually correct
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

            # Translate each element in node_list using mapping.map_name0
            translated_node_list = [mapping.label(item) for item in node_list]

            # Join the elements in node_list with "_"
            joined_node_string = "_".join(translated_node_list)

            # Add back the "COMPLEX:" prefix to the string
            complex_string = "COMPLEX:" + joined_node_string
        elif mapping.label(node):
            genesymbol = mapping.label(node)
            uniprot = mapping.id_from_label0(genesymbol)
        else:
            print("Error during translation, check syntax for ", node)

        if complex_string:
            new_entry = {"Genesymbol": complex_string, "Uniprot": node, "Type": "NaN"}
            self.nodes.loc[len(self.nodes)] = new_entry
        else:
            new_entry = {"Genesymbol": genesymbol, "Uniprot": uniprot, "Type": "NaN"}
            self.nodes.loc[len(self.nodes)] = new_entry

        self.nodes = self.nodes.drop_duplicates()
        return

    def add_edge(self, edge: pd.DataFrame):
        """
        Add an interaction to the list of interactions while converting it to the Omniflow-network format
        """
        # Check if the edge represents inhibition or stimulation and set the effect accordingly
        if (edge["is_inhibition"].values[0] == True and edge["is_stimulation"].values[0] == False):
            effect = "inhibition"
        elif (edge["is_inhibition"].values[0] == False and edge["is_stimulation"].values[0] == True):
            effect = "stimulation"
        else:
            effect = "undefined"

        # Get the type value from the edge DataFrame or set it to None
        edge_type = edge["type"].values[0] if "type" in edge.columns else None

        # Create a new DataFrame with edge information, including handling None for type
        df_edge = pd.DataFrame({
            "source": edge["source"],
            "target": edge["target"],
            "Type": edge_type,
            "Effect": effect
        })

        # add the new nodes to the nodes dataframe
        if not edge["source"].values[0] in self.nodes["Uniprot"].unique():
            self.add_node(edge["source"].values[0])
        if not edge["target"].values[0] in self.nodes["Uniprot"].unique():
            self.add_node(edge["target"].values[0])

        # Concatenate the new edge DataFrame with the existing edges in the graph
        self.edges = pd.concat([self.edges, df_edge])

    def add_paths_to_edge_list(self, paths):
        database = self.resources

        for path in paths:  # Iterate through the list of paths
            if isinstance(path, (str, tuple)):  # Handle single string or tuple
                path = [path]

            for i in range(0, len(path)):
                if i == len(path):
                    return
                interaction = database.loc[(database["source"] == path[i]) &
                                           (database["target"] == path[i - 1])]
                if not interaction.empty:
                    self.add_edge(interaction)

        return

    def convert_series_to_dataframes(self):
        return

    def connect_nodes(self):
        """
        Basic node connections. It adds all the interactions found in the omnipath database.
        Once an interaction is found it will be added to the list of edges
        """
        all_interactions = self.resources

        if len(self.nodes) == 1:
            print("Number of node insufficient to create connection")
        else:
            for node1, node2 in combinations(self.nodes["Uniprot"], 2):
                interaction1 = all_interactions.loc[(all_interactions["source"] == node1) &
                                                    (all_interactions["target"] == node2)]
                interaction1 = interaction1.reset_index(drop=True)
                interaction2 = all_interactions.loc[(all_interactions["source"] == node2) &
                                                    (all_interactions["target"] == node1)]
                interaction2 = interaction2.reset_index(drop=True)
                if not interaction1.empty:
                    self.add_edge(interaction1)
                if not interaction2.empty:
                    self.add_edge(interaction2)
        return

    def is_connected(self) -> bool:
        """
        Return True if all the nodes in the nodes list are connected, otherwise it returns False
        """
        return False

    def complete_connection(self,
                            maxlen: int = 2,
                            mode: Literal['OUT', 'IN', 'ALL'] = 'ALL'):
        """
        This function tries to connect all nodes of a network object using one of the methods presented in the Connection
        object (in the methods folder, enrichment_methods.py). This should be a core characteristic of this package and
        the user should have the possibility to choose different methods to enrich its Network object.

        TO COMPLETE STILL WORK IN PROGRESS
        """
        connect = Connections(self.resources)  # here we have the database where we look for the interactions

        connect_network = Connections(self.edges)  # here we have the edges dataframe of the network
        for node1, node2 in combinations(self.nodes["Uniprot"], 2):
            i = 0
            while i <= maxlen:
                print("looking for paths in the network with length: ", i, " for node ", node1, " and ", node2)
                paths = connect_network.find_paths(node1, node2, maxlen=i, mode="ALL")  # as first step, I make sure
                # that there is at least one path between two nodes in the network
                if not paths and i < maxlen:  # if there is no path, look for another one until reach maxlen
                    i += 1
                    continue
                elif not paths and i == maxlen:  # if there is no path of maxlen, look for a new path in the database
                    # until we find a new one
                    flag = False
                    i = 0
                    while not flag:
                        print("Looking for paths with length: ", i, " for node ", node1, " and ", node2)
                        paths = connect.find_paths(node1, node2, maxlen=i, mode=mode)
                        if not paths:
                            i += 1
                        else:
                            self.add_paths_to_edge_list(paths)
                            flag = True
                    break
                elif paths:  # if there is a path, there is no need to connect the node, so we iterate through another
                    # pair of nodes
                    break
        self.edges = self.edges.drop_duplicates()
        return

    def convert_edgelist_into_genesymbol(self):
        """
        This function converts the edge dataframe from uniprot to genesymbol. This may be removed later,
        I am not sure if it can be useful, mainly because the edge list will contain many different entities that
        will not be possible to translate, and so it will be useless...
        """
        # Convert the 'source' column of the edge DataFrame from UniProt to gene symbols
        self.edges["source"] = self.edges["source"].apply(lambda x: mapping.label(x) or x)

        # Convert the 'target' column of the edge DataFrame from UniProt to gene symbols
        self.edges["target"] = self.edges["target"].apply(lambda x: mapping.label(x) or x)

        return

