def fill(network,
                  only_signed: bool = False,
                  consensus_only: bool = False):
    """
    Basic node connections. It adds all the interactions found in the omnipath database.
    Once an interaction is found it will be added to the list of edges.
    The only_signed flag makes sure that just signed interaction will be added to the network, while "consensus_only"
    makes sure that just signed interaction with consensus among references will be included.

    Parameters:
    - only_signed: A boolean flag indicating whether to only add signed interactions to the network.
    - consensus_only: A boolean flag indicating whether to only add signed interactions with consensus among references to the network.

    Returns:
    None. The function modifies the network object in-place.
    """

    if len(network.nodes) == 1:
        print("Number of node insufficient to create connection")
        return

    def add_edge_if_not_empty_and_signed(node1, node2):
        """
        Helper function to add an edge to the network if the interaction is not empty and, if the `only_signed` flag is set, the interaction is signed.

        Parameters:
        - node1: The source node of the interaction.
        - node2: The target node of the interaction.

        Returns:
        None. The function modifies the network object in-place.
        """
        interaction = network.resources.loc[(network.resources["source"] == node1) &
                                         (network.resources["target"] == node2)]
        if not interaction.empty and (not only_signed or check_sign(interaction, consensus_only) != "undefined"):
            network.add_edge(interaction)

    for node1, node2 in combinations(network.nodes["Uniprot"], 2):
        add_edge_if_not_empty_and_signed(node1, node2)
        add_edge_if_not_empty_and_signed(node2, node1)
    return


