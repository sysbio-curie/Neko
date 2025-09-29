from typing import Any, Dict, List, Optional

import pandas as pd


class NetworkState:
    """Immutable snapshot of a network, with parent/child bookkeeping."""

    def __init__(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
        metadata: Optional[Dict[str, Any]] = None,
        state_id: Optional[int] = None,
        parent_ids: Optional[List[int]] = None,
    ):
        self.nodes = nodes.copy(deep=True)
        self.edges = edges.copy(deep=True)
        self.metadata: Dict[str, Any] = metadata or {}
        self.state_id = state_id
        self.parent_ids: List[int] = list(parent_ids or [])
        self.children_ids: List[int] = []

    def add_child(self, child_id: int) -> None:
        if child_id not in self.children_ids:
            self.children_ids.append(child_id)

    def __repr__(self) -> str:  # pragma: no cover - debug helper
        return (
            f"<NetworkState id={self.state_id} parents={self.parent_ids} "
            f"children={self.children_ids} metadata={self.metadata}>"
        )
