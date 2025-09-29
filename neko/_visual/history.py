from __future__ import annotations

from typing import TYPE_CHECKING, Iterable

import networkx as nx
from graphviz import Digraph

if TYPE_CHECKING:  # pragma: no cover
    from ..core.network import Network


def build_history_graph(network: "Network") -> nx.DiGraph:
    """Return a directed graph of the network state history."""

    graph = nx.DiGraph()
    for state_id in network.list_states():
        # list_states returns dictionaries with id/metadata
        graph.add_node(state_id["id"], **state_id.get("metadata", {}))
    for state_id, state in network._states.items():  # pylint: disable=protected-access
        for child in state.children_ids:
            graph.add_edge(state_id, child)
    return graph


def history_digraph(network: "Network", include_metadata: bool = True) -> Digraph:
    """Render history as a Graphviz Digraph."""

    digraph = Digraph()

    for state_id, state in network._states.items():  # pylint: disable=protected-access
        metadata = state.metadata or {}
        label_lines = [f"State {state_id}"]
        if metadata.get("label"):
            label_lines.append(str(metadata["label"]))
        if include_metadata and metadata.get("method"):
            args = metadata.get("args", [])
            kwargs = metadata.get("kwargs", {})
            params = ", ".join(args)
            if kwargs:
                kwargs_repr = ", ".join(f"{k}={v}" for k, v in kwargs.items())
                params = ", ".join(filter(None, [params, kwargs_repr]))
            label_lines.append(f"{metadata['method']}({params})")
        elif include_metadata:
            extras = {
                key: value for key, value in metadata.items()
                if key not in {"label", "method", "args", "kwargs"}
            }
            for key, value in extras.items():
                label_lines.append(f"{key}: {value}")
        digraph.node(str(state_id), "\n".join(label_lines))
        for child in state.children_ids:
            digraph.edge(str(state_id), str(child))
    return digraph


def history_html(network: "Network", include_metadata: bool = True, div_class: str = "neko-history-graph") -> str:
    """Return an HTML snippet embedding the history graph as SVG."""

    svg = history_digraph(network, include_metadata=include_metadata).pipe(format="svg").decode("utf-8")
    return f'<div class="{div_class}">{svg}</div>'
