"""
graph_traversal.py

This module contains graph traversal algorithms (DFS, BFS) for use in the NeKo package.
"""
from typing import List, Callable, Any

def dfs_algorithm(find_paths_func: Callable, node1: str, node2: str, maxlen: int, only_signed: bool, consensus: bool, connect_with_bias: bool, add_paths_func: Callable, connect_nodes_func: Callable, edges_df: Any) -> None:
    """
    Depth-First Search (DFS) algorithm to find paths between two nodes in the network.
    """
    paths = find_paths_func(start=node1, end=node2, maxlen=maxlen, minlen=1, only_signed=only_signed, consensus=consensus)
    if paths:
        add_paths_func(paths)
        if connect_with_bias:
            connect_nodes_func(only_signed, consensus)
            edges_df.drop_duplicates(inplace=True)
            edges_df.reset_index(drop=True, inplace=True)


def bfs_algorithm(bfs_func: Callable, node1: str, node2: str, maxlen: int, only_signed: bool, consensus: bool, connect_with_bias: bool, add_paths_func: Callable, connect_nodes_func: Callable, edges_df: Any) -> None:
    """
    Breadth-First Search (BFS) algorithm to find paths between two nodes in the network.
    """
    paths = bfs_func(start=node1, end=node2, maxlen=maxlen, only_signed=only_signed, consensus=consensus)
    if paths:
        add_paths_func(paths)
        if connect_with_bias:
            connect_nodes_func(only_signed, consensus)
            edges_df.drop_duplicates(inplace=True)
            edges_df.reset_index(drop=True, inplace=True)

