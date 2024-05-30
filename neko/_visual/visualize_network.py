from graphviz import Digraph
from IPython.display import display


def wrap_node_name(node_name):
    if ":" in node_name:
        node_name = node_name.replace(":", "_")
    if node_name.startswith("COMPLEX"):
        return node_name[8:]
    else:
        return node_name


class NetworkVisualizer:
    def __init__(self, network, predefined_node=None, color_by="Effect"):
        net = network.copy()
        net.convert_edgelist_into_genesymbol()
        self.dataframe_edges = net.edges
        self.dataframe_nodes = net.nodes
        self.color_by = color_by
        self.graph = Digraph(format='pdf')
        self.edge_colors = {
            'stimulation': 'green',
            'inhibition': 'red',
            'form complex': 'blue',
            'undefined': 'gray',  # Adding the mapping for "undefined" effect
            'bimodal': 'darkorange',  # Adding the mapping for "bimodal" effect
            # Add more custom mappings if needed
        }
        self.node_colors = {}  # Dictionary to store custom node colors
        self.predefined_node = wrap_node_name(predefined_node) if predefined_node else None
        # Apply wrap_node_name function to node names in dataframe_nodes
        self.dataframe_nodes['Genesymbol'] = self.dataframe_nodes['Genesymbol'].apply(wrap_node_name)
        self.dataframe_nodes['Uniprot'] = self.dataframe_nodes['Uniprot'].apply(wrap_node_name)

        # Apply wrap_node_name function to node names in "source" and "target" columns of dataframe_edges
        self.dataframe_edges['source'] = self.dataframe_edges['source'].apply(wrap_node_name)
        self.dataframe_edges['target'] = self.dataframe_edges['target'].apply(wrap_node_name)

    def set_custom_edge_colors(self, custom_edge_colors):
        # Update the edge_colors dictionary with custom mappings
        self.edge_colors.update(custom_edge_colors)

    def set_node_colors(self, node_colors):
        # Update the node_colors dictionary with custom node colorsdataframe_nodes
        self.node_colors.update(node_colors)

    def add_edges_to_graph(self):
        for _, row in self.dataframe_edges.iterrows():
            effect = row['Effect']
            source = wrap_node_name(row['source'])
            target = wrap_node_name(row['target'])

            # Display only edges connected to the predefined node
            if self.predefined_node and (source != self.predefined_node and target != self.predefined_node):
                continue

            # Determine edge attributes based on effect
            if effect == 'stimulation':
                arrowhead = 'normal'
                color = 'green'
                dir = 'forward'
            elif effect == 'inhibition':
                arrowhead = 'tee'
                color = 'red'
                dir = 'forward'
            elif effect == 'form complex':
                arrowhead = 'dot'
                color = 'blue'
                dir = 'forward'
            elif effect == 'bimodal':
                arrowhead = 'box'
                color = 'darkorange'
                dir = 'forward'
            else:
                arrowhead = 'normal'  # Default arrow
                color = 'black'
                dir = 'none'

            # Add the edge to the graph with specified attributes
            self.graph.edge(source, target, color=color, arrowhead=arrowhead, dir=dir)

    def add_nodes_to_graph(self):
        for _, row in self.dataframe_nodes.iterrows():
            node = row['Genesymbol']
            # add function to set color
            node_color = self.node_colors.get(node, 'lightgray')

            # Display only the predefined node and its connections
            if self.predefined_node and (node != self.predefined_node):
                continue

            wrapped_node = wrap_node_name(node)
            self.graph.node(wrapped_node, style='filled', fillcolor=node_color)

    def tissue_mapping(self, tissue_df):
        """
        Color the nodes based on their expression in the tissue of interest (based on data from The Human Protein Atlas).

        Args:
            tissue_df (DataFrame): DataFrame containing results indicating whether each gene symbol has tissue annotations containing the selected tissue.
        """
        for _, row in tissue_df.iterrows():
            gene_symbol = row['Genesymbol']
            in_tissue = row['in_tissue']
            node_color = 'lightblue' if in_tissue else 'lightgray'
            self.node_colors[gene_symbol] = node_color

    def render(self, output_file='network', view=False):
        self.add_edges_to_graph()
        self.add_nodes_to_graph()

        if view:
            self.graph.view(filename=output_file)
        else:
            self.graph.render(filename=output_file)
            display(self.graph)  # Display the graph directly in the Jupyter Notebook
