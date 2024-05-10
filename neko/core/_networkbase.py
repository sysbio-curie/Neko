import pandas as pd


class NetworkBase:


    def __init__(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
    ):

        self.nodes = nodes
        self.edges = edges
