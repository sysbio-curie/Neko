from graphviz import Digraph
import pandas as pd


class NetworkVisualizer:
    def __init__(self, dataframe, color_by="Effect"):
        self.dataframe = dataframe
        self.color_by = color_by
        self.graph = Digraph(format='pdf')

    def add_edges_to_graph(self):
        for _, row in self.dataframe.iterrows():
            effect = row['Effect']
            edge_color = 'green' if effect == 'stimulation' else 'red'

            if self.color_by == "type":
                edge_color = row['Type']

            self.graph.edge(row['source'], row['target'], color=edge_color,
                            dir='forward' if effect == 'stimulation' else 'none')

    def render(self, output_file='network', view=False):
        self.add_edges_to_graph()

        if view:
            self.graph.view(filename=output_file)
        else:
            self.graph.render(filename=output_file)
