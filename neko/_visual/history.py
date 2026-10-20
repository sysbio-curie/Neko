from __future__ import annotations

from typing import TYPE_CHECKING, Iterable

import networkx as nx
from graphviz import Digraph
import html
import json
from ..core.history import json_value

if TYPE_CHECKING:  # pragma: no cover
    from ..core.network import Network


def build_history_graph(network: "Network") -> nx.DiGraph:
    """Return a directed graph of the network state history."""

    graph = nx.DiGraph()
    graph.graph.update(database_name=network.database_name,
                       current_state_id=network.current_state_id,
                       pruned_state_ids=list(network._pruned_state_ids))
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
        if state_id == network.current_state_id:
            label_lines.append("current")
        if include_metadata and metadata.get("method"):
            args = metadata.get("args", [])
            kwargs = metadata.get("kwargs", {})
            params = ", ".join(args)
            if kwargs:
                kwargs_repr = ", ".join(f"{k}={v}" for k, v in kwargs.items())
                params = ", ".join(filter(None, [params, kwargs_repr]))
            label_lines.append(f"{metadata['method']}({params[:160]})")
        elif include_metadata:
            extras = {
                key: value for key, value in metadata.items()
                if key in {"description", "database_name"}
            }
            for key, value in extras.items():
                label_lines.append(f"{key}: {value}")
        if metadata.get('ancestry_gaps'):
            label_lines.append('ancestry gap (pruned states)')
        digraph.node(str(state_id), "\n".join(label_lines),
                     URL=f"#neko-state-{state_id}" if include_metadata else "",
                     tooltip=f"State {state_id}: structured explanation")
        for child in state.children_ids:
            digraph.edge(str(state_id), str(child))
    return digraph


def history_html(network: "Network", include_metadata: bool = True, div_class: str = "neko-history-graph") -> str:
    """Return an HTML snippet embedding the history graph as SVG."""

    svg = history_digraph(network, include_metadata=include_metadata).pipe(format="svg").decode("utf-8")
    details = ''
    if include_metadata:
        details = ''.join(
            f'<details id="neko-state-{sid}"><summary>State {sid} explanation</summary>'
            f'<pre>{html.escape(json.dumps(json_value(state.metadata), indent=2))}</pre></details>'
            for sid, state in network._states.items()
        )
    return f'<div class="{html.escape(div_class, quote=True)}">{svg}{details}</div>'
