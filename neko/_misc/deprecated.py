
def check_gene_list_format(gene_list: list[str]) -> list[str]:
    """
    This function checks the format of the gene list and returns True if the gene list is in Uniprot format, False if the gene list is in genesymbol format.

    Parameters:
    - gene_list: A list of gene identifiers. The gene identifiers can be either Uniprot identifiers or genesymbols.

    Returns:
    - A boolean indicating whether the gene list is in Uniprot format (True) or genesymbol format (False).
    """
    # Check if the gene list contains Uniprot identifiers
    if all(mapping.id_from_label0(gene) for gene in gene_list):
        return True
    # Check if the gene list contains genesymbols
    elif all(mapping.label(gene) for gene in gene_list):
        return False

