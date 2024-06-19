def phenotype(network,
                               phenotype: str = None,
                               id_accession: str = None,
                               sub_genes: list[str] = None,
                               maxlen: int = 2,
                               only_signed: bool = False,
                               compress: bool = False
                               ):
    """
    This function connects genes to a phenotype based on the provided parameters. It retrieves phenotype markers,
    identifies unique Uniprot genes, and connects them to the network. It also has the option to compress the network
    by substituting specified genes with the phenotype name.

    Parameters:
    - phenotype: The phenotype to connect to. If not provided, it will be retrieved using the id_accession.
    - id_accession: The accession id of the phenotype. If not provided, the phenotype parameter must be given.
    - sub_genes: A list of genes to be considered for connection. If not provided, all nodes in the network are considered.
    - maxlen: The maximum length of the paths to be searched for.
    - only_signed: A boolean flag to indicate whether to filter unsigned paths.
    - compress: A boolean flag to indicate whether to substitute the specified genes with the phenotype name.

    Returns:
    None. The function modifies the network object in-place.
    """
    # Initialize lists for Uniprot and genesymbol genes
    uniprot_gene_list = []
    genesymbols_genes = []

    # Retrieve phenotype markers
    network._ensure_go()
    phenotype_genes = network.ontology.get_markers(phenotype=phenotype, id_accession=id_accession)
    if not phenotype_genes:
        print("Something went wrong while getting the markers for: ", phenotype, " and ", id_accession)
        print("Check URL and try again")
        return
    # Convert phenotype genes to Uniprot identifiers
    uniprot_genes = [translate_id(i)[2] for i in phenotype_genes]
    # If sub_genes are provided, check their format and convert to Uniprot or genesymbol as needed
    if sub_genes:
        if check_gene_list_format(sub_genes):
            uniprot_gene_list = sub_genes
            genesymbols_genes = [translate_id(i)[2] for i in sub_genes]
        else:
            uniprot_gene_list = [translate_id(i)[2] for i in sub_genes]
            genesymbols_genes = sub_genes

    print("Starting connecting network's nodes to: ", phenotype_genes)
    # Identify unique Uniprot genes not already in the network
    unique_uniprot = set(uniprot_genes) - set(uniprot_gene_list if uniprot_gene_list else network.nodes["Uniprot"])
    # Identify unique genesymbols not already in the network
    unique_genesymbol = set(phenotype_genes) - set(
        genesymbols_genes if genesymbols_genes else network.nodes["Genesymbol"])
    # Connect the network's nodes to the unique Uniprot genes
    network.connect_component(uniprot_gene_list if uniprot_gene_list else network.nodes["Uniprot"].tolist(),
                           list(unique_uniprot),
                           mode="OUT",
                           maxlen=maxlen, only_signed=only_signed)

    # If compress is True, substitute specified genes with the phenotype name
    if compress:
        phenotype = phenotype or network.ontology.accession_to_phenotype_dict[id_accession]
        phenotype_modified = phenotype.replace(" ", "_")

        # Substitute the specified genes with the phenotype name in the nodes dataframe
        network.nodes['Uniprot'] = network.nodes['Uniprot'].apply(
            lambda x: phenotype_modified if x in unique_uniprot else x)
        network.nodes['Genesymbol'] = network.nodes['Genesymbol'].apply(
            lambda x: phenotype_modified if x in unique_genesymbol else x)

        # Substitute the specified genes with the phenotype name in the edges dataframe
        for column in ['source', 'target']:
            network.edges[column] = network.edges[column].apply(
                lambda x: phenotype_modified if x in unique_uniprot else x)

        # Group by source and target, and aggregate with the custom function for each column
        network.edges = network.edges.groupby(['source', 'target']).agg({
            'Type': join_unique,  # Aggregate types with the custom function
            'Effect': join_unique,  # Aggregate effects with the custom function
            'References': join_unique  # Aggregate references with the custom function
        }).reset_index()

        # Identify common genes between Uniprot genes and the network's nodes
        common_genes = set(uniprot_genes).intersection(
            set(uniprot_gene_list if uniprot_gene_list else network.nodes["Uniprot"]))
        # For each common gene, add a new edge connecting the gene to the phenotype
        for gene in common_genes:
            new_edge = {"source": gene, "target": phenotype_modified, "Effect": "stimulation",
                        "References": "Gene Ontology"}
            network.edges = network.edges.append(new_edge, ignore_index=True)
        network.edges.drop_duplicates()
    return

