def filter_unsigned_paths(network, paths: list[tuple], consensus: bool) -> list[tuple]:
    """
    This function filters out unsigned paths from the provided list of paths. An unsigned path is a path where at least
    one interaction does not have a defined sign (stimulation or inhibition). The function checks each interaction in each
    path, and if the interaction is unsigned, the path is not included in the output list.

    Parameters:
    - paths: A list of paths, where each path is a sequence of nodes.
    - consensus: A boolean indicating whether to check for consensus among references when determining the sign of an interaction.

    Returns:
    - A list of paths where all interactions in each path are signed.
    """
    interactions = network.resources
    filtered_paths = []
    for path in paths:
        is_full_signed = True
        for i in range(0, len(path)):
            if i == len(path) - 1:
                break
            interaction = interactions.loc[(interactions["source"] == path[i]) &
                                           (interactions["target"] == path[i + 1])]
            if not interaction.empty and check_sign(interaction, consensus) == "undefined":
                is_full_signed = False
        if is_full_signed:
            filtered_paths.append(path)
    return filtered_paths


def check_sign(
        interaction: pd.DataFrame | dict,
        consensus: bool = False
    ) -> Literal['stimulation', 'inhibition', 'form_complex', 'undefined']:
    """
    This function checks the sign of an interaction in the Omnipath format (Pandas DataFrame or Series).
    The attribute "consensus" checks for the consistency of the sign of the interaction among the references.

    Parameters:
    - interaction: A pandas DataFrame or Series representing the interaction.
    - consensus: A boolean indicating whether to check for consensus among references.

    Returns:
    - A string indicating the sign of the interaction: "stimulation", "inhibition", "form complex", or "undefined".
    """
    # Handle both DataFrame and Series input
    if isinstance(interaction, pd.DataFrame):
        interaction = interaction.iloc[0]

    signs = ('stimulation', 'inhibition', 'form_complex', 'undefined')
    prefix = 'consensus' if consensus else 'is'

    for sign in signs:

        if interaction.get(f'{prefix}_{sign}', False):

            return sign


def add_cascade_to_edge_list(network, cascades):
    """
    This function adds cascades to the edge list of the network. A cascade is a sequence of nodes where each node is
    connected to the next node in the sequence. The function checks if there is an interaction between each pair of nodes
    in the cascade in the resources' database. If an interaction exists, it is added to the edge list of the network.

    Parameters:
    - cascades: A list of cascades, where each cascade is a sequence of nodes.

    Returns:
    None. The function modifies the network object in-place.
    """
    database = network.resources

    for cascade in cascades:
        interaction_in = database.loc[(database["source"] == cascade[0]) &
                                      (database["target"] == cascade[1])]
        if interaction_in.empty:
            print("Empty interaction for node ", cascade[0], " and ", cascade[1])
        else:
            network.add_edge(interaction_in)
    network.edges = network.edges.drop_duplicates()

    return


def add_paths_to_edge_list(network, paths):
    """
    This method adds paths to the edge list of the network. A path is a sequence of nodes where each node is
    connected to the next node in the sequence. The function checks if there is an interaction between each pair of nodes
    in the path in the resources' database. If an interaction exists, it is added to the edge list of the network.

    Parameters:
    - paths: A list of paths, where each path is a sequence of nodes. A node can be a string or a tuple.

    Returns:
    None. The function modifies the network object in-place by adding the interactions to the edges DataFrame.
    """
    # Access the resources database
    database = network.resources

    # Iterate through the list of paths
    for path in paths:
        # Handle single string or tuple
        if isinstance(path, (str, tuple)):
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
                network.add_edge(interaction)

    # Remove duplicate edges from the edge list
    network.edges = network.edges.drop_duplicates()

    return


def remove_path(network, path: list[str]):
    """
    This function removes a path from the network. It takes a list of nodes representing the path and removes all the edges
    between the nodes in the path.

    Parameters:
    - path: A list of nodes representing the path to be removed. The nodes can be represented as strings or tuples.

    Returns:
    None. The function modifies the network object in-place by removing the edges between the nodes in the path.
    """
    # check if node1 and node2 are in genesymbol format or uniprot format
    if check_gene_list_format(path):
        path = [translate_id(node)[2] for node in path]

    # Iterate through the nodes in the path and remove the edges between them
    for i in range(0, len(path)):
        if i == len(path) - 1:
            break
        network.remove_edge(path[i], path[i + 1])
    return


def drop_missing_nodes(network):
    """
    This function drops the nodes that are not present in the resources database and print a warning with the name of the missing nodes.

    The function works as follows:
    1. It first calls the `check_nodes` function to get a list of nodes that exist in the resources' database.
    2. It then finds the nodes in the network that are not in this list, and removes them from the network.
    3. If there are any missing nodes, it prints a warning with their names.

    This function does not return anything. It modifies the `nodes` attribute of the `Network` object in-place.
    """
    # Get the list of nodes that exist in the resources database
    existing_nodes = network.check_nodes(network.nodes["Uniprot"].tolist())

    # Find the nodes in the network that are not in the list of existing nodes
    missing_nodes = [node for node in network.nodes["Uniprot"].tolist() if node not in existing_nodes]

    # Remove the missing nodes from the network
    network.nodes = network.nodes[~network.nodes["Uniprot"].isin(missing_nodes)]

    # Print a warning with the name of the missing nodes
    if missing_nodes:
        print(
            "Warning: The following nodes were not found in the resources database and have been removed from the "
            "network:",
            ", ".join(missing_nodes))
    return
