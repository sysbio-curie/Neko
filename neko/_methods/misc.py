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
