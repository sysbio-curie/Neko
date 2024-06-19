

def complete_connection(network,
                        maxlen: int = 2,
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

    # Set the search depth based on the k_mean parameter
    if k_mean == 'tight':
        i_search = 3
    elif k_mean == 'extensive':
        i_search = 4

    # Create a Connections object for the resources
    connect = Connections(network.resources)

    # Copy the nodes
    nodes = network.nodes.copy()

    # Create a Connections object for the edges
    connect_network = Connections(network.edges)

    # Iterate through all combinations of nodes
    for node1, node2 in tqdm(combinations(nodes["Uniprot"], 2), desc="Connecting nodes"):
        print("Connecting nodes %s and %s" % (node1, node2))
        if not network.check_node_existence(node1) or not network.check_node_existence(node2):
            print(
                "Error: node %s is not present in the resources database" % node1 if not network.check_node_existence(
                    node1) else "Error: node %s is not present in the resources database" % node2)
            continue
        i = 0
        # Reset the object connect_network, updating the possible list of paths if minimal is True
        if minimal:
            connect_network = Connections(network.edges)

        # Start a loop to find paths between nodes
        while i <= maxlen:
            # As first step, make sure that there is at least one path between two nodes in the network
            paths_in = connect_network.find_paths(node2, node1, maxlen=i)
            paths_out = connect_network.find_paths(node1, node2, maxlen=i)

            # If paths in both directions are found, break the loop
            if paths_in and paths_out:
                break

            # If no paths are found and the maximum length has not been reached, increment the length and continue the loop
            if not (paths_in and paths_out) and i < maxlen:
                i += 1
                continue

            # If no paths are found and the maximum length has been reached, search for new paths in the database
            if not paths_in and i == maxlen:
                flag = False
                j = 0
                # Start a loop to find new paths in the database
                while not flag and j <= i_search:
                    print("Searching for paths from %s to %s with length %s" % (node2, node1, j))
                    paths_in = connect.find_paths(node2, node1, maxlen=j)
                    # Filter unsigned paths if only_signed is True
                    if only_signed:
                        paths_in = network.filter_unsigned_paths(paths_in, consensus)
                    # If no paths are found, increment the length and continue the loop
                    if not paths_in:
                        j += 1
                    else:
                        # If a path is found, add it to the edge list and connect nodes if connect_node_when_first_introduced is True
                        network.add_paths_to_edge_list(paths_in)
                        if connect_node_when_first_introduced:
                            network.connect_nodes(only_signed, consensus)
                            network.edges = network.edges.drop_duplicates()
                        flag = True

            # Repeat the same process for paths in the opposite direction
            if not paths_out and i == maxlen:
                flag = False
                j = 0
                while not flag and j <= i_search:
                    print("Searching for paths from %s to %s with length %s" % (node1, node2, j))
                    paths_out = connect.find_paths(node1, node2, maxlen=j)
                    if only_signed:
                        paths_out = network.filter_unsigned_paths(paths_out, consensus)
                    if not paths_out:
                        j += 1
                    else:
                        network.add_paths_to_edge_list(paths_out)
                        if connect_node_when_first_introduced:
                            network.connect_nodes(only_signed, consensus)
                            network.edges = network.edges.drop_duplicates()
                        flag = True
            break

    # If connect_node_when_first_introduced is False, connect nodes after all paths have been found
    if not connect_node_when_first_introduced:
        network.connect_nodes(only_signed, consensus)
        network.edges = network.edges.drop_duplicates()
    return

def complete_connection_2(network,
                          maxlen: int = 2,
                          k_mean: int = None,
                          minimal: bool = True,
                          only_signed: bool = False,
                          consensus: bool = False,
                          connect_node_when_first_introduced: bool = True):
    """
    This function attempts to connect all nodes of a network object using one of the methods presented in the Connection
    object. This is a core characteristic of this package and the user should have the possibility to choose different
    methods to enrich its Network object.

    Parameters:
    - maxlen: The maximum length of the paths to be searched for. Default is 2.
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

    # Create a Connections object for the resources
    connect = Connections(network.resources)

    # Copy the nodes
    nodes = network.nodes.copy()

    # Create a Connections object for the edges
    connect_network = Connections(network.edges)

    i = 0
    while i <= maxlen:
        # Iterate through all combinations of nodes
        for node1, node2 in tqdm(combinations(nodes["Uniprot"], 2), desc="Connecting nodes"):
            print("Connecting nodes %s and %s" % (node1, node2))
            if not network.check_node_existence(node1) or not network.check_node_existence(node2):
                print(
                    "Error: node %s is not present in the resources database" % node1 if not network.check_node_existence(
                        node1) else "Error: node %s is not present in the resources database" % node2)
                continue

            # Reset the object connect_network, updating the possible list of paths if minimal is True
            if minimal:
                connect_network = Connections(network.edges)

            paths_in = []
            paths_out = []

            # As first step, make sure that there is at least one path between two nodes in the network within the k mean
            if k_mean:
                for j in range(k_mean):
                    paths_in = connect_network.find_paths(node2, node1, maxlen=j)
                    paths_out = connect_network.find_paths(node1, node2, maxlen=j)
            else:
                paths_in = connect_network.find_paths(node2, node1, maxlen=i)
                paths_out = connect_network.find_paths(node1, node2, maxlen=i)

            # If paths in both directions are found, break the loop
            if paths_in and paths_out:
                continue

            # If no paths are found and the maximum length has been reached, search for new paths in the database
            if not paths_out:
                print("Searching for paths from %s to %s with length %s" % (node1, node2, i))
                paths_out = connect.find_paths(node1, node2, maxlen=i)
                # Filter unsigned paths if only_signed is True
                if only_signed:
                    paths_out = network.filter_unsigned_paths(paths_out, consensus)
                # If no paths are found, increment the length and continue the loop
                if paths_out:
                    # If a path is found, add it to the edge list and connect nodes if connect_node_when_first_introduced is True
                    network.add_paths_to_edge_list(paths_out)
                    if connect_node_when_first_introduced:
                        network.connect_nodes(only_signed, consensus)
                        network.edges = network.edges.drop_duplicates()

            # If no paths are found and the maximum length has been reached, search for new paths in the database
            if not paths_in:
                print("Searching for paths from %s to %s with length %s" % (node2, node1, i))
                paths_in = connect.find_paths(node2, node1, maxlen=i)
                # Filter unsigned paths if only_signed is True
                if only_signed:
                    paths_in = network.filter_unsigned_paths(paths_in, consensus)
                # If no paths are found, increment the length and continue the loop
                if paths_in:
                    # If a path is found, add it to the edge list and connect nodes if connect_node_when_first_introduced is True
                    network.add_paths_to_edge_list(paths_in)
                    if connect_node_when_first_introduced:
                        network.connect_nodes(only_signed, consensus)
                        network.edges = network.edges.drop_duplicates()

        i = i + 1

    # If connect_node_when_first_introduced is False, connect nodes after all paths have been found
    if not connect_node_when_first_introduced:
        network.connect_nodes(only_signed, consensus)
        network.edges = network.edges.drop_duplicates()
    return

