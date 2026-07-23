"""
Connection strategies for NeKo networks.

This module contains high-level strategies for connecting nodes in a Network object.
Each function should accept a Network instance as the first argument.
"""
import logging
from itertools import combinations
from typing import List, Optional, Union

import pandas as pd
from typing_extensions import Literal

from .._methods.enrichment_methods import Connections
from .tools import consolidate_edges, is_connected

logger = logging.getLogger(__name__)


def _node_identifiers(nodes: pd.DataFrame) -> pd.Series:
    """Return each node's edge identifier, with a label fallback."""

    return nodes["Uniprot"].where(
        nodes["Uniprot"].notna(),
        nodes["Genesymbol"],
    )


def _remove_nodes(network, nodes: pd.DataFrame) -> bool:
    """Remove rows and report whether the node count decreased."""

    previous_count = len(network.nodes)

    for node in _node_identifiers(nodes).dropna().tolist():
        network.remove_node(node)

    made_progress = len(network.nodes) < previous_count

    if not made_progress:
        logger.warning(
            'Stopped disconnected-node cleanup because no removable '
            'identifier was found.',
        )

    return made_progress


def connect_nodes(network, only_signed: bool = False, consensus_only: bool = False) -> None:
    """
    Basic node connections. Adds all interactions found in the resources database.
    """
    if len(network.nodes) == 1:
        print("Number of node insufficient to create connection")
        return

    def add_edge_if_not_empty_and_signed(node1, node2):
        if node2 in network._connect.find_all_neighbours(node1):
            interaction = network.resources.loc[(network.resources["source"] == node1) &
                                                (network.resources["target"] == node2)]
            if not interaction.empty and (
                not only_signed or network.check_sign(interaction, consensus_only) != "undefined"):
                network.add_edge(interaction)

    for node1, node2 in combinations(network.nodes["Uniprot"], 2):
        add_edge_if_not_empty_and_signed(node1, node2)
        add_edge_if_not_empty_and_signed(node2, node1)
    return

def connect_subgroup(network, group, maxlen: int = 1, only_signed: bool = False, consensus: bool = False) -> None:
    """
    Connect all nodes in a subgroup by finding paths between all pairs and adding them to the network.
    """
    if not network.check_gene_list_format(group):
        uniprot_gene_list = group
    else:
        uniprot_gene_list = [network.mapping_node_identifier(i)[2] for i in group]
    if len(uniprot_gene_list) == 1:
        print("Number of node insufficient to create connection")
    else:
        for node1, node2 in combinations(uniprot_gene_list, 2):
            i = 0
            paths_in = []
            paths_out = []
            while i <= maxlen:
                if not paths_out:
                    paths_out = network._connect.find_paths(node1, node2, maxlen=i, only_signed=only_signed, consensus=consensus)
                if not paths_in:
                    paths_in = network._connect.find_paths(node2, node1, maxlen=i, only_signed=only_signed, consensus=consensus)
                if (not paths_in or not paths_out) and i <= maxlen:
                    i += 1
                if ((paths_in or paths_out) and i > maxlen) or (paths_in and paths_out):
                    paths = paths_out + paths_in
                    network._add_paths_to_edge_list(paths)
                    break
    return

def connect_component(network, comp_A, comp_B, maxlen: int = 2, mode: Literal['OUT', 'IN', 'ALL'] = 'OUT', only_signed: bool = False, consensus: bool = False) -> None:
    """
    Connect subcomponents of a network using the specified mode and add paths to the network.
    """
    if mode == "IN":
        paths = network._connect.find_paths(comp_B, comp_A, maxlen=maxlen, only_signed=only_signed, consensus=consensus)
    elif mode == "OUT":
        paths = network._connect.find_paths(comp_A, comp_B, maxlen=maxlen, only_signed=only_signed, consensus=consensus)
    elif mode == "ALL":
        paths = network._connect.find_paths(comp_A, comp_B, maxlen=maxlen, only_signed=only_signed, consensus=consensus) + \
                network._connect.find_paths(comp_B, comp_A, maxlen=maxlen, only_signed=only_signed, consensus=consensus)
    else:
        print("The only accepted modes are IN, OUT or ALL, please check the syntax")
        return
    network._add_paths_to_edge_list(paths)
    all_nodes = set(network.nodes['Uniprot'].values)
    set_a = set(comp_A)
    set_b = set(comp_B)
    set_c = list(all_nodes.difference(set_a).difference(set_b))
    if len(set_c) > 0:
        connect_subgroup(network, set_c, only_signed=only_signed, maxlen=maxlen, consensus=consensus)
    return

def connect_to_upstream_nodes(network, nodes_to_connect=None, depth: int = 1, rank: int = 1, only_signed: bool = True, consensus: bool = False) -> None:
    """
    Connect provided nodes to their upstream nodes in the network.
    """
    if nodes_to_connect is None:
        nodes_to_connect = network.nodes["Uniprot"].tolist()
    cascades = network._connect.find_upstream_cascades(nodes_to_connect, depth, rank)
    # No need to filter cascades, as sign filtering is handled in Connections if needed
    network._add_cascade_to_edge_list(cascades)
    network.edges.drop_duplicates().reset_index(drop=True)
    return

def connect_genes_to_phenotype(
        network,
        phenotype: str = None,
        id_accession: str = None,
        sub_genes: list = None,
        maxlen: int = 2,
        only_signed: bool = False,
        compress: bool = False,
        taxon_id=9606,
        include_descendants: bool = False,
        exclude_automatic_assertions: bool = False,
    ) -> None:
    """
    Connect a network to GO-associated genes and optionally compress them.

    GO-provided UniProt identifiers are used directly. Gene-symbol mapping is
    retained only as a fallback for associations in another identifier space.
    """
    id_accession = network._ontology.resolve_accession(
        phenotype=phenotype,
        id_accession=id_accession,
    )
    term = network._ontology.get_term(id_accession)
    go_genes = network._ontology.fetch_go_genes(
        id_accession,
        taxon_id=taxon_id,
        include_descendants=include_descendants,
        exclude_automatic_assertions=exclude_automatic_assertions,
    )
    if not go_genes:
        logger.warning(
            "No genes associated with %s for taxon %s.",
            term.go_id,
            taxon_id,
        )
        return

    uniprot_genes = []
    for gene in go_genes:
        uniprot = None
        if gene.gene_id and gene.gene_id.startswith("UniProtKB:"):
            uniprot = gene.gene_id.split(":", 1)[1]
        if uniprot is None and gene.symbol:
            uniprot = network.mapping_node_identifier(gene.symbol)[2]
        if uniprot is None:
            logger.warning(
                "Skipping GO gene without a usable network identifier: %s",
                gene.gene_id or gene.symbol,
            )
            continue
        uniprot_genes.append(uniprot)

    if not uniprot_genes:
        logger.warning(
            "No genes associated with %s could be mapped to network IDs.",
            term.go_id,
        )
        return

    uniprot_gene_list = []
    if sub_genes:
        for gene in sub_genes:
            _, _, uniprot = network.mapping_node_identifier(gene)
            uniprot_gene_list.append(uniprot or gene)

    source_uniprot = (
        uniprot_gene_list
        if uniprot_gene_list
        else network.nodes["Uniprot"].dropna().tolist()
    )
    unique_uniprot = set(uniprot_genes) - set(source_uniprot)
    connect_component(
        network,
        source_uniprot,
        sorted(unique_uniprot),
        mode="OUT",
        maxlen=maxlen,
        only_signed=only_signed,
    )
    if compress:
        phenotype_modified = term.label.replace(" ", "_")
        compressed_node_mask = network.nodes['Uniprot'].isin(unique_uniprot)
        has_compressed_nodes = compressed_node_mask.any()

        network.nodes = network.nodes.loc[~compressed_node_mask].copy()
        for column in ['source', 'target']:
            network.edges[column] = network.edges[column].replace(
                dict.fromkeys(unique_uniprot, phenotype_modified),
            )

        common_genes = set(uniprot_genes).intersection(source_uniprot)
        if has_compressed_nodes or common_genes:
            phenotype_node = pd.DataFrame([{
                'Genesymbol': phenotype_modified,
                'Uniprot': phenotype_modified,
                'Type': 'phenotype',
            }])
            network.nodes = pd.concat(
                [network.nodes, phenotype_node],
                ignore_index=True,
            )
            network.nodes = network.nodes.drop_duplicates(
                subset=['Genesymbol', 'Uniprot'],
            ).reset_index(drop=True)

        if common_genes:
            go_edges = pd.DataFrame({
                'source': sorted(common_genes),
                'target': phenotype_modified,
                'Type': 'gene ontology association',
                'Effect': 'stimulation',
                'References': f'Gene Ontology: {term.go_id}',
            })
            network.edges = pd.concat(
                [network.edges, go_edges],
                ignore_index=True,
            )

        network.edges = consolidate_edges(network.edges)

        if hasattr(network, 'sync_nodes_from_df'):
            network.sync_nodes_from_df()
        if hasattr(network, 'sync_edges_from_df'):
            network.sync_edges_from_df()
    return

def connect_network_radially(network, max_len: int = 1, direction: Literal['OUT', 'IN', None] = None, loops: bool = False, consensus: bool = False, only_signed: bool = True) -> None:
    """
    Connect all nodes of a network in a radial manner.
    """
    initial_nodes = network.initial_nodes
    initial_nodes_set = {
        identifier
        for node in initial_nodes
        if (identifier := network.mapping_node_identifier(node)[2]) is not None
    }
    i = 0
    source_nodes = initial_nodes_set
    target_nodes = initial_nodes_set
    while i < max_len:
        new_nodes = []
        if direction == 'OUT' or direction is None:
            for source in source_nodes:
                target_neighs = network._connect.find_target_neighbours(source)
                if source in target_neighs and not loops:
                    target_neighs.remove(source)
                target_paths = [(source, node) for node in target_neighs if (not only_signed or network._connect.is_signed_edge(source, node, consensus))]
                network._add_paths_to_edge_list(target_paths)
                target_neighs_filtered = [path[1] for path in target_paths]
                target_neighs_filtered = [node for node in target_neighs_filtered if node not in initial_nodes_set]
                new_nodes.extend(target_neighs_filtered)
            source_nodes = new_nodes
        new_nodes = []
        if direction == 'IN' or direction is None:
            for target in target_nodes:
                source_neighs = network._connect.find_source_neighbours(target)
                if target in source_neighs and not loops:
                    source_neighs.remove(target)
                source_paths = [(node, target) for node in source_neighs if (not only_signed or network._connect.is_signed_edge(node, target, consensus))]
                network._add_paths_to_edge_list(source_paths)
                source_neighs_filtered = [path[0] for path in source_paths]
                source_neighs_filtered = [node for node in source_neighs_filtered if node not in initial_nodes_set]
                new_nodes.extend(source_neighs_filtered)
            target_nodes = new_nodes
        i += 1
    # Remove disconnected nodes
    target_nodes_set = set(network.edges["target"].unique())
    source_nodes_set = set(network.edges["source"].unique())
    node_identifiers = _node_identifiers(network.nodes)
    disconnected_nodes = network.nodes[
        ~node_identifiers.isin(initial_nodes_set) & (
            ~node_identifiers.isin(target_nodes_set)
            | ~node_identifiers.isin(source_nodes_set)
        )
    ]
    while not disconnected_nodes.empty:
        if not _remove_nodes(network, disconnected_nodes):
            break
        target_nodes_set = set(network.edges["target"].unique())
        source_nodes_set = set(network.edges["source"].unique())
        node_identifiers = _node_identifiers(network.nodes)
        disconnected_nodes = network.nodes[
            ~node_identifiers.isin(initial_nodes_set) & (
                ~node_identifiers.isin(target_nodes_set)
                | ~node_identifiers.isin(source_nodes_set)
            )
        ]
    return

def connect_as_atopo(network, strategy: Literal['radial', 'complete', None] = None, max_len: int = 1, loops: bool = False, outputs=None, only_signed: bool = True, consensus: bool = False) -> None:
    """
    Connect all nodes of a network in a topological manner.
    """
    initial_nodes = [
        identifier
        for node in network.initial_nodes
        if (identifier := network.mapping_node_identifier(node)[2]) is not None
    ]
    initial_nodes_set = set(initial_nodes)
    if strategy == 'radial':
        connect_network_radially(network, max_len, direction=None, loops=loops, consensus=consensus, only_signed=only_signed)
    elif strategy == 'complete':
        network.complete_connection(max_len, minimal=True, only_signed=only_signed, consensus=consensus, connect_with_bias=False)
    starting_nodes = set(_node_identifiers(network.nodes).dropna())
    if outputs is None:
        return

    outputs_uniprot = []
    invalid_outputs = []

    for node in outputs:
        if network.add_node(node):
            identifier = network.mapping_node_identifier(node)[2]

            if identifier is not None:
                outputs_uniprot.append(identifier)
                continue

        invalid_outputs.append(node)

    if invalid_outputs:
        logger.warning(
            'Ignoring output nodes without a usable resource identifier: %s',
            ', '.join(map(str, invalid_outputs)),
        )

    if not outputs_uniprot:
        logger.warning('No valid output nodes were available for connection.')
        return

    depth = 1
    while not is_connected(network):
        connect_to_upstream_nodes(network, outputs_uniprot, depth=depth, rank=len(outputs_uniprot), only_signed=only_signed, consensus=consensus)
        new_nodes = set(network.nodes["Uniprot"].tolist()) - starting_nodes
        new_nodes = new_nodes - set(outputs_uniprot)
        for node in new_nodes:
            if node not in network.edges["target"].unique():
                network.remove_node(node)
            if not loops and any((network.edges['source'] == node) & (network.edges['target'] == node)):
                network.remove_node(node)
        if depth == 4:
            print("Current depth is 4, stopping the process")
            break
        depth += 1
    network.edges.drop_duplicates().reset_index(drop=True)
    target_nodes_set = set(network.edges["target"].unique())
    node_identifiers = _node_identifiers(network.nodes)
    disconnected_nodes = network.nodes[
        ~node_identifiers.isin(initial_nodes_set)
        & ~node_identifiers.isin(target_nodes_set)
    ]
    while not disconnected_nodes.empty:
        if not _remove_nodes(network, disconnected_nodes):
            break
        target_nodes_set = set(network.edges["target"].unique())
        node_identifiers = _node_identifiers(network.nodes)
        disconnected_nodes = network.nodes[
            ~node_identifiers.isin(initial_nodes_set)
            & ~node_identifiers.isin(target_nodes_set)
        ]
    return

def complete_connection(network,
                        maxlen: Optional[int] = 2,
                        algorithm: Literal['bfs', 'dfs'] = 'dfs',
                        minimal: bool = True,
                        only_signed: bool = False,
                        consensus: bool = False,
                        connect_with_bias: bool = False,
                        ) -> None:
    """
    Attempts to connect all nodes of a network object using one of the methods presented in the Connection object.
    For each node pair, checks for existing paths in both directions using BFS (with sign/consensus as needed).
    If a path is missing, calls the selected algorithm to try to find and add a path.
    Uses the minimal flag to reset the Connections object as needed.
    After all, if connect_with_bias is False, calls connect_nodes and deduplicates edges.
    If maxlen=None, performs a single unbounded BFS (no iterative deepening needed).
    """
    nodes = network.nodes.copy()
    connect_network = Connections(network.edges)

    for node1, node2 in combinations(nodes["Uniprot"], 2):
        if not network.check_node(node1) or not network.check_node(node2):
            continue
        if minimal:
            connect_network = Connections(network.edges)
        def find_bfs_path(start, end):
            return connect_network.bfs(start=start, end=end, maxlen=maxlen, only_signed=only_signed, consensus=consensus)
        # Check for existing paths in both directions
        paths_in = find_bfs_path(node2, node1)
        paths_out = find_bfs_path(node1, node2)
        # If a path is missing, call the selected algorithm to try to find and add a path
        if not paths_in:
            network._algorithms[algorithm](node1=node2, node2=node1, maxlen=maxlen, only_signed=only_signed, consensus=consensus, connect_with_bias=connect_with_bias)
            connect_network = Connections(network.edges)
        if not paths_out:
            network._algorithms[algorithm](node1=node1, node2=node2, maxlen=maxlen, only_signed=only_signed, consensus=consensus, connect_with_bias=connect_with_bias)
            connect_network = Connections(network.edges)
    if not connect_with_bias:
        network.connect_nodes(only_signed, consensus)
        network.edges = network.edges.drop_duplicates().reset_index(drop=True)
    return
