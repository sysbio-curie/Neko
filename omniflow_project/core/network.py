from pypath.utils import mapping
import omnipath as op
from itertools import combinations
import pandas as pd
from .._inputs.resources import Resources
from .._methods.enrichment_methods import Connections


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
        self.edges = pd.DataFrame(columns=["Source", "Target", "Type", "Effect"])
        self.resources = Resources().all_omnipath_interactions()
        if initial_nodes:
            for node in initial_nodes:
                self.add_node(node)

    def add_node(self,
                 node: str):
        """
        Add a node to the list of nodes checking that the syntax for the genesymbol is actually correct
        """
        uniprot = mapping.map_name0(node, 'genesymbol', 'uniprot')
        genesymbol = mapping.label(uniprot)
        new_entry = {"Genesymbol": genesymbol, "Uniprot": uniprot, "Type": "NaN"}
        self.nodes.loc[len(self.nodes)] = new_entry

        return

    def add_edge(self,
                 edge: pd.DataFrame):

        if (edge["is_inhibition"].values[0] == True and edge["is_stimulation"].values[0] == False):
            effect = "inhibition"
        elif (edge["is_inhibition"].values[0] == False and edge["is_stimulation"].values[0] == True):
            effect = "stimulation"
        else:
            effect = "undefined"
        df_edge = pd.DataFrame({"Source": edge["source"], "Target": edge["target"], "Type": edge["type"],
                                "Effect": effect})
        self.edges = pd.concat([self.edges, df_edge])

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

    def complete_connection(self):
        """
        This function tries to connect all nodes of a network object using one of the methods presented in the Connection
        object (in the methods folder, enrichment_methods.py). This should be a core characteristic of this package and
        the user should have the possibility to choose different methods to enrich its Network object.

        TO COMPLETE STILL WORK IN PROGRESS
        """
        connect = Connections()
        all_interactions = self.resources
        edges = Connections.force_nodes_connections(self.nodes, self.edges, self.resources)
        return


