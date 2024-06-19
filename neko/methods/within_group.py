def within_group(network,
                     group: (str | pd.DataFrame | list[str]),
                     maxlen: int = 1,
                     only_signed: bool = False,
                     consensus: bool = False
                     ):
    """
    This function is used to connect all the nodes in a particular subgroup. It iterates over all pairs of nodes in the
    subgroup and finds paths between them in the resources' database. If a path is found, it is added to the edge list of
    the network. The function also filters out unsigned paths if the `only_signed` flag is set to True.

    Parameters:
    - group: A list of nodes representing the subgroup to connect. Nodes can be represented as strings, pandas DataFrame, or list of strings.
    - maxlen: The maximum length of the paths to be searched for in the resources' database. Default is 1.
    - only_signed: A boolean flag indicating whether to only add signed interactions to the network. Default is False.
    - consensus: A boolean flag indicating whether to only add signed interactions with consensus among references to the network. Default is False.

    Returns:
    None. The function modifies the network object in-place.
    """
    connect = Connections(network.resources)  # here we have the database where we look for the interactions
    if not check_gene_list_format(group):
        uniprot_gene_list = group
    else:
        uniprot_gene_list = [translate_id(i)[2] for i in group]
    if len(uniprot_gene_list) == 1:
        print("Number of node insufficient to create connection")
    else:
        for node1, node2 in combinations(uniprot_gene_list, 2):
            i = 0
            paths_in = []
            paths_out = []
            while i <= maxlen:
                if not paths_out:
                    paths_out = connect.find_paths(node1, node2, maxlen=i)
                    if only_signed:
                        paths_out = network.filter_unsigned_paths(paths_out, consensus)
                if not paths_in:
                    paths_in = connect.find_paths(node2, node1, maxlen=i)
                    if only_signed:
                        paths_in = network.filter_unsigned_paths(paths_in, consensus)
                if not paths_in or not paths_out and i <= maxlen:
                    i += 1
                if (paths_in or paths_out) and i > maxlen or (paths_in and paths_out):
                    paths = paths_out + paths_in
                    network.add_paths_to_edge_list(paths)
                    break
    return


