from __future__ import annotations

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


def check_sign(interaction: pd.DataFrame,
               consensus: bool = False) -> str:
    """
    This function quickly check the sign of an interaction in the omnipath format (Pandas Series)
    The attribute "consensus" check for the consistency of the sign of the interaction among the references.
    """
    if consensus:
        if interaction["consensus_stimulation"].values[0]:
            return "stimulation"
        elif interaction["consensus_inhibition"].values[0]:
            return "inhibition"
        else:
            return "undefined"
    else:
        if interaction["is_stimulation"].values[0]:
            return "stimulation"
        elif interaction["is_inhibition"].values[0]:
            return "inhibition"
        elif interaction["form_complex"].values[0]:
            return "form complex"
        else:
            return "undefined"


def mapping_node_identifier(node: str) -> list[str]:
    """
    Returns an array with the possible identifier of the node
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

    return [complex_string, genesymbol, uniprot]


def translate_paths(paths):
    translated_list = []

    # If input_list is a list of strings
    if isinstance(paths[0], str):
        translated_list = [mapping_node_identifier(item)[1] or mapping_node_identifier(item)[0] for item
                           in paths]
    # If input_list is a list of lists of strings
    elif isinstance(paths[0], list):
        for sublist in paths:
            translated_sublist = [mapping_node_identifier(item)[1] or mapping_node_identifier(item)[0] for
                                  item in sublist]
            translated_list.append(translated_sublist)

    return translated_list


def join_unique(series):
    # Filter out None values before converting to set and joining
    filtered_series = [str(item) for item in series if item is not None]
    unique_items = set(filtered_series)
    return ', '.join(unique_items)


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
            print("Loading deafault omnipath all interactions")
            res = Resources()
            res.load_omnipath_interactions()
            self.resources = res.interactions  ### in future release, a string can determine which database to use/load
        if initial_nodes:
            for node in initial_nodes:
                self.add_node(node)
        elif sif_file:
            self.load_network_from_sif(sif_file)

    def copy(self):
        new_instance = copy.deepcopy(self)
        return new_instance

    def add_node(self,
                 node: str):
        """
        Add a node to the list of nodes checking that the syntax for the genesymbol is actually correct
        """
        complex_string, genesymbol, uniprot = mapping_node_identifier(node)

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
        print(edge["source"].values[0], edge["target"].values[0])
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
        """Load a network object from sif file"""
        interactions = []
        node_set = set()

        with open(sif_file, "r") as f:
            for line in f:
                if line.startswith('#'):  # Skip comment lines
                    continue
                interaction = line.strip().split()
                if len(interaction) < 3:
                    continue  # Skip malformed lines

                # Determine the effect based on the interaction type
                if interaction[1] in ["1", "activate", "stimulate", "phosphorilate", "stimulation"]:
                    effect = "stimulation"
                elif interaction[1] in ["-1", "inhibit", "block", "inhibition"]:
                    effect = "inhibition"
                elif interaction[1] in ["form complex", "form_complex", "form-complex", "complex formation"]:
                    effect = "form complex"
                else:
                    effect = "undefined"

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
        database = self.resources

        for path in paths:  # Iterate through the list of paths
            if isinstance(path, (str, tuple)):  # Handle single string or tuple
                path = [path]

            for i in range(0, len(path)):
                if i == len(path) - 1:
                    break
                interaction_in = database.loc[(database["source"] == path[i]) &
                                              (database["target"] == path[i + 1])]
                interaction_out = database.loc[(database["target"] == path[i]) &
                                               (database["source"] == path[i + 1])]

                if not interaction_in.empty and (mode == "ALL" or "IN"):
                    self.add_edge(interaction_in)
                elif not interaction_out.empty and (mode == "ALL" or "OUT"):
                    self.add_edge(interaction_out)
        self.edges = self.edges.drop_duplicates()
        return

    def convert_series_to_dataframes(self):
        return

    def connect_nodes(self,
                      only_signed: bool = False,
                      consensus_only: bool = False):
        """
        Basic node connections. It adds all the interactions found in the omnipath database.
        Once an interaction is found it will be added to the list of edges.
        The only_signed flag makes sure that just signed interaction will be added to the network, while "consensus_only"
        makes sure that just signed interaction with consensus among references will be included.
        """
        all_interactions = self.resources

        if len(self.nodes) == 1:
            print("Number of node insufficient to create connection")
        else:
            for node1, node2 in combinations(self.nodes["Uniprot"], 2):
                interaction1 = all_interactions.loc[(all_interactions["source"] == node1) &
                                                    (all_interactions["target"] == node2)]
                # interaction1 = interaction1.reset_index(drop=True)
                interaction2 = all_interactions.loc[(all_interactions["source"] == node2) &
                                                    (all_interactions["target"] == node1)]
                # interaction2 = interaction2.reset_index(drop=True)
                if only_signed:
                    if not interaction1.empty and check_sign(interaction1, consensus_only) != "undefined":
                        self.add_edge(interaction1)
                    if not interaction2.empty and check_sign(interaction2, consensus_only) != "undefined":
                        self.add_edge(interaction2)
                else:
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

    def filter_unsigned_paths(self,
                              paths: list[tuple],
                              consensus: bool) -> list[tuple]:
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
                         group: (
                             str | pd.DataFrame | list[str]
                         ),
                         maxlen: int = 1,
                         mode: Literal['OUT', 'IN', 'ALL'] = 'ALL',
                         only_signed: bool = False,
                         consensus: bool = False
                         ):
        """
        Function used to connect all the nodes in a particular subgroup.
        """

        connect = Connections(self.resources)  # here we have the database where we look for the interactions
        if len(self.nodes) == 1:
            print("Number of node insufficient to create connection")
        else:
            for node1, node2 in combinations(group, 2):
                flag = False
                i = 0
                while not flag and i <= maxlen:
                    print("looking for paths in the database with length: ", i, " for node ",
                          mapping_node_identifier(node1)[1], " and ", mapping_node_identifier(node2)[1])
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
                        print(translate_paths(paths))
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
        This function tries to connect all nodes of a network object using one of the methods presented in the Connection
        object (in the methods folder, enrichment_methods.py). This should be a core characteristic of this package and
        the user should have the possibility to choose different methods to enrich its Network object.

        TO COMPLETE STILL WORK IN PROGRESS
        """
        if k_mean == 'tight':
            i_search = 3
        elif k_mean == 'extensive':
            i_search = 4

        connect = Connections(self.resources)  # here we have the database where we look for the interactions
        nodes = self.nodes.copy()
        connect_network = Connections(self.edges)  # here we have the edges dataframe of the network
        for node1, node2 in combinations(nodes["Uniprot"], 2):
            i = 0
            if minimal:
                connect_network = Connections(self.edges)  # here we have the edges dataframe of the network
                # basically with this move I am re-setting the object connect_network, updating the possible list of
                # paths. In this way, next time I am checking if two nodes are already connected in my network,
                # I am including also the newly found genes/protein! this can reduce by a lot the final net size
            while i <= maxlen:
                print("looking for paths in the network with length: ", i, " for node ",
                      mapping_node_identifier(node1)[1], " and ", mapping_node_identifier(node2)[1])
                # as first step, I make sure that there is at least one path between two nodes in the network
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
                    print(translate_paths(paths))
                if not paths and i < maxlen:  # if there is no path, look for another one until reach maxlen
                    i += 1
                    continue
                elif not paths and i == maxlen:  # if there is no path of maxlen, look for a new path in the database
                    # until we find a new one
                    flag = False
                    i = 0
                    while not flag and i <= i_search:
                        print("i_search = ", i_search)
                        print("Looking for paths with length: ", i, " for node ",
                              mapping_node_identifier(node1)[1], " and ", mapping_node_identifier(node2)[1])
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
                            print(translate_paths(paths))
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
                              str | list[str]
                          ),
                          comp_B: (
                              str | list[str]
                          ),
                          maxlen: int = 2,
                          mode: Literal['OUT', 'IN', 'ALL'] = 'OUT',
                          only_signed: bool = False,
                          consensus: bool = False
                          ):
        """
        This function tries to connect subcomponents of a network object using one of the methods presented in the Connection
        object (in the methods folder, enrichment_methods.py). This should be a core characteristic of this package and
        the user should have the possibility to choose different methods to enrich its Network object.

        TO COMPLETE STILL WORK IN PROGRESS
        """
        print(comp_A)
        print(comp_B)
        connect = Connections(self.resources)
        if mode == "IN":
            paths_in = connect.find_paths(comp_A, comp_B, maxlen=maxlen, mode=mode)
            paths = paths_in
        elif mode == "OUT":
            paths_out = connect.find_paths(comp_A, comp_B, maxlen=maxlen, mode=mode)
            paths = paths_out
        elif mode == "ALL":
            paths_out = connect.find_paths(comp_A, comp_B, maxlen=maxlen, mode="OUT")
            paths_in = connect.find_paths(comp_A, comp_B, maxlen=maxlen, mode="IN")
            paths = paths_out + paths_in
        else:
            print("The only accepted modes are IN, OUT or ALL, please check the syntax")
            return
        if only_signed:
            paths = self.filter_unsigned_paths(paths, consensus)
        self.add_paths_to_edge_list(paths, mode)
        all_nodes = set(self.nodes['Uniprot'].values)
        set_a = set(comp_A)
        set_b = set(comp_B)
        set_c = all_nodes.difference(set_a)
        set_c = set_c.difference(set_b)
        set_c = list(set_c)
        print(set_c)
        self.connect_subgroup(set_c, only_signed=only_signed, maxlen=maxlen)
        return

    def convert_edgelist_into_genesymbol(self):
        """
        This function converts the edge dataframe from uniprot to genesymbol. This may be removed later,
        I am not sure if it can be useful, mainly because the edge list will contain many different entities that
        will not be possible to translate, and so it will be useless...
        """
        self.edges.loc[:, "source"] = self.edges["source"].apply(
            lambda x: mapping_node_identifier(x)[0] or mapping_node_identifier(x)[1])

        self.edges.loc[:, "target"] = self.edges["target"].apply(
            lambda x: mapping_node_identifier(x)[0] or mapping_node_identifier(x)[1])

        return

    def connect_genes_to_phenotype(self,
                                   phenotype: str = None,
                                   id_accession: str = None,
                                   sub_genes: list[str] = None,
                                   maxlen: int = 2,
                                   only_signed: bool = False,
                                   compress: bool = False
                                   ):
        phenotype_genes = self.ontology.get_markers(phenotype=phenotype, id_accession=id_accession)
        if not phenotype_genes:
            print("Something went wrong while getting the markers for: ", phenotype, " and ", id_accession)
            print("Check URL and try again")
            return
        uniprot_genes = [mapping_node_identifier(i)[2] for i in phenotype_genes]

        print("Starting connecting network's nodes to: ", phenotype_genes)
        if sub_genes:
            unique_uniprot = set(uniprot_genes) - set(sub_genes)
            unique_uniprot_list = list(unique_uniprot)
            self.connect_component(sub_genes, unique_uniprot_list, mode="OUT", maxlen=maxlen, only_signed=only_signed)
        else:
            nodes = list(self.nodes["Uniprot"])
            unique_uniprot = set(uniprot_genes) - set(nodes)
            unique_uniprot_list = list(unique_uniprot)
            self.connect_component(nodes, unique_uniprot_list, mode="OUT", maxlen=maxlen, only_signed=only_signed)

        if compress:
            if not phenotype and id_accession:
                phenotype = self.ontology.accession_to_phenotype_dict[id_accession]

            phenotype_modified = phenotype.replace(" ", "_")

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
