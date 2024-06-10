import pandas as pd

def compare_networks(network1, network2):

    """
    This function compares two networks and identifies common, unique, and conflicting interactions and nodes.
    Args:
        network1:
        network2:

    Returns:

    """

    df1 = network1.convert_edgelist_into_genesymbol()
    df2 = network2.convert_edgelist_into_genesymbol()
    # Merge dataframes on 'source' and 'target' to identify common and conflicting interactions
    merged = pd.merge(df1, df2, on=['source', 'target'], suffixes=('_1', '_2'), how='outer', indicator=True)

    # Determine the comparison type for interactions
    merged['comparison'] = merged['_merge'].apply(
        lambda x: 'Unique to Network 1' if x == 'left_only' else ('Unique to Network 2' if x == 'right_only' else 'Common')
    )

    # Convert the comparison column to categorical and add the new category
    merged['comparison'] = pd.Categorical(merged['comparison'], categories=['Unique to Network 1', 'Unique to Network 2', 'Common', 'Conflicting'])

    # Update the comparison type for conflicting interactions
    conflicting_mask = (merged['_merge'] == 'both') & (merged['Effect_1'] != merged['Effect_2'])
    merged.loc[conflicting_mask, 'comparison'] = 'Conflicting'

    # Create the interaction comparison DataFrame
    interaction_comparison = merged[['source', 'target', 'comparison']]

    # Determine unique and common nodes
    nodes_1 = set(df1['source']).union(set(df1['target']))
    nodes_2 = set(df2['source']).union(set(df2['target']))

    unique_nodes_network_1 = nodes_1 - nodes_2
    unique_nodes_network_2 = nodes_2 - nodes_1
    common_nodes = nodes_1 & nodes_2

    # Create a list of node comparison data
    node_comparison_data = []
    for node in unique_nodes_network_1:
        node_comparison_data.append([node, 'Unique to Network 1'])
    for node in unique_nodes_network_2:
        node_comparison_data.append([node, 'Unique to Network 2'])
    for node in common_nodes:
        node_comparison_data.append([node, 'Common'])

    # Create the node comparison DataFrame
    node_comparison = pd.DataFrame(node_comparison_data, columns=['node', 'comparison'])

    return interaction_comparison, node_comparison
