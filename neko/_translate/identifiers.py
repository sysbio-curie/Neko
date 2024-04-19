def translate_id(node: str) -> list[str]:
    """
    This function takes a node identifier and returns a list containing the possible identifiers for the node.
    The identifiers include a complex string, a genesymbol, and a uniprot identifier. The function uses the
    mapping.id_from_label0 and mapping.label functions from the pypath.utils.mapping module to translate the node
    identifier into these different formats.

    Parameters:
    - node: A string representing the node identifier. The node identifier can be a genesymbol, a uniprot identifier,
            or a complex string.

    Returns:
    - A list containing the complex string, genesymbol, and uniprot identifier for the node. If the node identifier
      cannot be translated into one of these formats, the corresponding value in the list is None.
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

        # Translate each element in node_list using mapping.label
        translated_node_list = [mapping.label(mapping.id_from_label0(item)) for item in node_list]

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
    """
    This function translates a list of paths, where each path is a sequence of node identifiers.
    It uses the helper function `handle_complex_identifier` to translate each node identifier in the paths.

    Parameters:
    - paths: A list of paths, where each path is a sequence of node identifiers.
             A node identifier can be a string or a list of strings.

    Returns:
    - A list of translated paths, where each path is a sequence of translated node identifiers.
    """
    translated_list = []

    def handle_complex_identifier(item):
        """
        This helper function translates a node identifier using the `translate_id` function.
        It checks all possible identifiers (complex, genesymbol, uniprot) and returns the first non-None value.

        Parameters:
        - item: A node identifier.

        Returns:
        - The translated node identifier.
        """
        identifiers = translate_id(item)
        return identifiers[0] or identifiers[1] or identifiers[2]

    # If input_list is a list of strings
    if isinstance(paths[0], str):
        translated_list = [handle_complex_identifier(item) for item in paths]
    # If input_list is a list of lists of strings
    elif isinstance(paths[0], list):
        for sublist in paths:
            translated_sublist = [handle_complex_identifier(item) for item in sublist]
            translated_list.append(translated_sublist)

    return translated_list



def convert_edgelist_into_genesymbol(network):
    """
    This function converts the edge dataframe from uniprot to genesymbol.
    """

    def convert_identifier(x):
        identifiers = translate_id(x)
        return identifiers[0] or identifiers[1]

    network.edges["source"] = network.edges["source"].apply(convert_identifier)
    network.edges["target"] = network.edges["target"].apply(convert_identifier)

    return
