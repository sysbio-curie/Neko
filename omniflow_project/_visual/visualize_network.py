from graphviz import Digraph
from IPython.display import display

class NetworkVisualizer:
    def __init__(self, network, color_by="Effect"):
        self.dataframe_edges = network.edges
        self.dataframe_nodes = network.nodes
        self.color_by = color_by
        self.graph = Digraph(format='pdf')
        self.edge_colors = {
            'stimulation': 'green',
            'inhibition': 'red',
            'undefined': 'gray',  # Adding the mapping for "undefined" effect
            # Add more custom mappings if needed
        }
        self.node_colors = {}  # Dictionary to store custom node colors

    def set_custom_edge_colors(self, custom_edge_colors):
        # Update the edge_colors dictionary with custom mappings
        self.edge_colors.update(custom_edge_colors)

    def set_node_colors(self, node_colors):
        # Update the node_colors dictionary with custom node colorsdataframe_nodes
        self.node_colors.update(node_colors)

    def add_edges_to_graph(self):
        for _, row in self.dataframe_edges.iterrows():
            effect = row['Effect']
            edge_color = self.edge_colors.get(effect, 'black')

            if self.color_by == "type":
                edge_color = self.edge_colors.get(row['Type'], 'black')

            self.graph.edge(row['source'], row['target'], color=edge_color,
                            dir='forward' if effect == 'stimulation' else 'none')

    def add_nodes_to_graph(self):
        for _, row in self.dataframe_nodes.iterrows():
            node = row['Genesymbol']
            node_color = self.node_colors.get(node, 'lightgray')
            self.graph.node(node, style='filled', fillcolor=node_color)

    def render(self, output_file='network', view=False):
        self.add_edges_to_graph()
        self.add_nodes_to_graph()

        if view:
            self.graph.view(filename=output_file)
        else:
            self.graph.render(filename=output_file)
            display(self.graph)  # Display the graph directly in the Jupyter Notebook
