from graphviz import Digraph
from IPython.display import display
from yfiles_jupyter_graphs import GraphWidget

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
    
    def yfiles_visual(
     #       self,
            network,
            graph_layout,
            directed,
        ):
        # creating empty object for visualization
        w = GraphWidget()

        objects = []
        for idx, item in network.nodes.iterrows():
            obj = {
                "id" : network.nodes["Uniprot"].loc[idx],
                "properties" : {"label": network.nodes["Genesymbol"].loc[idx]},
                "color":"#ffffff",
                "styles":{"backgroundColor":"#ffffff"}
                }
            objects.append(obj)
        w.nodes=objects

        # filling w with edges
        objects= []
        for index, row in network.edges.iterrows():
            obj={
            "id":network.edges["Effect"].loc[index],
            "start" : network.edges["source"].loc[index],
            "end" : network.edges["target"].loc[index],
            "properties":{"references":network.edges["References"].loc[index]}}
            objects.append(obj)
        w.edges=objects

        def custom_edge_color_mapping(edge: 'Dict'):
            """let the edge be purple if the starting node has an even index"""
            return ("#ff0066" if edge['id']  == "inhibition" else "#00cc00")
        w.set_edge_color_mapping(custom_edge_color_mapping)


        def custom_label_styles_mapping(node: 'Dict'):
            """let the label be the negated purple big index"""
            return {
                'backgroundColor': '#ffffff', 
                'color': '#ffffff', 
                'shape':'round-rectangle'
            }
        w.set_node_styles_mapping(custom_label_styles_mapping)

        w.directed=directed
        w.graph_layout=graph_layout

        w.show()

    def vis_comparison(
            int_comparison,
            node_comparison,
            graph_layout,
            directed,
        ):
        # creating empty object for visualization
        w = GraphWidget()

        objects = []
        for idx, item in node_comparison.iterrows():
            obj = {
                "id" : node_comparison["node"].loc[idx],
                "properties" : {"label": node_comparison["node"].loc[idx],
                               "comparison": node_comparison["comparison"].loc[idx],},
           #     "color":"#ffffff",
         #       "styles":{"backgroundColor":"#ffffff"}
                }
            objects.append(obj)
        w.nodes=objects

        # filling w with edges
        objects= []
        for index, row in int_comparison.iterrows():
            obj={
            "id":int_comparison["comparison"].loc[index],
            "start" : int_comparison["source"].loc[index],
            "end" : int_comparison["target"].loc[index]
                }
            objects.append(obj)
        w.edges=objects
    
        def custom_node_color_mapping(node: 'Dict'):
            if node['properties']['comparison'] == "Unique to Network 1":
                return "#ff0066"
            elif node['properties']['comparison'] == "Unique to Network 2":
                return "#00cc00"
            elif node['properties']['comparison'] == "Common":
                return "#0000ff"
        w.set_node_styles_mapping(custom_node_color_mapping)
    
        def custom_edge_color_mapping(edge: 'Dict'):
            if edge['id'] == "Unique to Network 1":
                return "#ff0066"
            elif edge['id'] == "Unique to Network 2":
                return "#00cc00"
            elif edge['id'] == "Common":
                return "#0000ff"
            elif edge['id'] == "Conflicting":
                return "#ffcc00"    
        w.set_edge_color_mapping(custom_edge_color_mapping)


        def custom_label_styles_mapping(node: 'Dict'):
            return {
                'backgroundColor': '#ffffff', 
                'color': '#ffffff', 
                'shape':'round-rectangle'
            }
        w.set_node_styles_mapping(custom_label_styles_mapping)
        w.directed=directed
        w.graph_layout=graph_layout

        w.show()    


