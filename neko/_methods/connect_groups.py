def connect_groups(network,
                      comp_A: (str | list[str]),
                      comp_B: (str | list[str]),
                      maxlen: int = 2,
                      mode: Literal['OUT', 'IN', 'ALL'] = 'OUT',
                      only_signed: bool = False,
                      consensus: bool = False):
    """
    This function attempts to connect subcomponents of a network object using one of the methods presented in the Connection
    object. This is a core characteristic of this package and the user should have the possibility to choose different
    methods to enrich its Network object.

    Parameters:
    - comp_A: A string or list of strings representing the first component to connect.
    - comp_B: A string or list of strings representing the second component to connect.
    - maxlen: The maximum length of the paths to be searched for. Default is 2.
    - mode: The search mode, which can be 'OUT', 'IN', or 'ALL'. Default is 'OUT'.
    - only_signed: A boolean flag indicating whether to filter unsigned paths. Default is False.
    - consensus: A boolean flag indicating whether to check for consensus among references. Default is False.

    Returns:
    None. The function modifies the network object in-place.

    """
    # Print the components for debugging purposes
    print(comp_A)
    print(comp_B)

    # Create a Connections object for the resources
    connect = Connections(network.resources)

    # Determine the search mode and find paths accordingly
    if mode == "IN":
        paths_in = connect.find_paths(comp_B, comp_A, maxlen=maxlen)
        paths = paths_in
    elif mode == "OUT":
        paths_out = connect.find_paths(comp_A, comp_B, maxlen=maxlen)
        paths = paths_out
    elif mode == "ALL":
        paths_out = connect.find_paths(comp_A, comp_B, maxlen=maxlen)
        paths_in = connect.find_paths(comp_B, comp_A, maxlen=maxlen)
        paths = paths_out + paths_in
    else:
        print("The only accepted modes are IN, OUT or ALL, please check the syntax")
        return

    # Filter unsigned paths if the only_signed flag is set
    if only_signed:
        paths = network.filter_unsigned_paths(paths, consensus)

    # Add the paths to the edge list
    network.add_paths_to_edge_list(paths)

    # Create sets of nodes for each component and the entire network
    all_nodes = set(network.nodes['Uniprot'].values)
    set_a = set(comp_A)
    set_b = set(comp_B)

    # Find the nodes that are not in either component
    set_c = all_nodes.difference(set_a)
    set_c = set_c.difference(set_b)
    set_c = list(set_c)

    # If there are nodes not in either component, connect them as a subgroup
    if len(set_c) > 0:
        network.connect_subgroup(set_c, only_signed=only_signed, maxlen=maxlen)

    return
