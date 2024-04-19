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

