from typing import Any, Dict, List, Optional
import pandas as pd

class NetworkState:
    """
    Represents a snapshot of a Network's state (nodes, edges, and metadata) at a given point in time.
    Intended for supporting undo/redo, branching, and provenance tracking in the Network class.

    Args:
        nodes (pd.DataFrame): The nodes DataFrame at this state.
        edges (pd.DataFrame): The edges DataFrame at this state.
        metadata (Optional[Dict[str, Any]]): Optional metadata (e.g., timestamp, description, provenance info).
    """
    def __init__(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        self.nodes = nodes.copy(deep=True)
        self.edges = edges.copy(deep=True)
        self.metadata = metadata or {}

    def __repr__(self):
        return f"<NetworkState nodes={len(self.nodes)}, edges={len(self.edges)}, metadata={self.metadata}>"
