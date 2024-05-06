def common_regulators(network,
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
        connect = Connections(network.resources)
        if nodes_to_connect is None:
            nodes_to_connect = network.nodes["Uniprot"].tolist()

        cascades = connect.find_upstream_cascades(nodes_to_connect, depth, rank)
        if only_signed:
            cascades = network.filter_unsigned_paths(cascades, consensus)
        network.add_cascade_to_edge_list(cascades)
        network.edges.drop_duplicates()
    except Exception as e:
        print(f"An error occurred while connecting to upstream nodes: {e}")
    return


