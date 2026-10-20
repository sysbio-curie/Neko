"""Structured records for the existing NetworkState graph.

These helpers do not maintain a second history. JSON records are deliberately
independent of provider data: restoring a snapshot is not construction replay.
"""

import copy
import inspect
import json
from functools import wraps
from pathlib import Path

import pandas as pd

from .strategy_options import UNSET


def json_value(value):
    """Convert supported scientific values without display truncation."""
    if value is UNSET:
        return {"unspecified": True}
    if isinstance(value, pd.DataFrame):
        return frame_payload(value)
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, set):
        return sorted((json_value(v) for v in value), key=canonical)
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if hasattr(value, "item"):
        return json_value(value.item())
    if isinstance(value, float):
        return None if pd.isna(value) else value
    if value is pd.NA or value is pd.NaT:
        return None
    raise TypeError(f"Unsupported history value: {type(value).__name__}")


def canonical(value):
    return json.dumps(json_value(value), sort_keys=True, allow_nan=False)


def copy_frame(frame):
    # DataFrame.copy(deep=True) does not copy nested dictionaries/lists.
    result = frame.copy(deep=True)
    for column in result.select_dtypes(include="object"):
        for index, value in result[column].items():
            if isinstance(value, (dict, list, set, tuple)):
                result.at[index, column] = copy.deepcopy(value)
    return result


def records(frame):
    return json_value(frame.to_dict("records"))


def frame_payload(frame):
    return {"columns": list(frame.columns), "records": records(frame)}


def frame_from_payload(payload):
    if (
        not isinstance(payload, dict)
        or not isinstance(payload.get("columns"), list)
        or not all(isinstance(column, str) for column in payload["columns"])
        or len(set(payload["columns"])) != len(payload["columns"])
        or not isinstance(payload.get("records"), list)
        or any(
            not isinstance(row, dict) or set(row) != set(payload["columns"])
            for row in payload["records"]
        )
    ):
        raise ValueError(
            "Malformed history table; columns and records must agree."
        )
    return pd.DataFrame(
        copy.deepcopy(payload["records"]), columns=payload["columns"]
    )


def frame_fingerprint(frame):
    return canonical(frame_payload(frame))


def unique_edges(frame):
    if frame.empty:
        return frame.copy().reset_index(drop=True)
    keys = pd.Series([canonical(row) for row in records(frame)])
    return frame.iloc[(~keys.duplicated()).to_numpy()].reset_index(drop=True)


def provenance(row, origin="unknown", database_name=None):
    existing = row.get("Provenance")
    if isinstance(existing, dict):
        return copy.deepcopy(existing)
    return {
        "origin": origin,
        "original": {
            "origin": origin,
            "database_name": database_name if origin == "resource" else None,
            "effect": json_value(row.get("Effect")),
            "references": json_value(row.get("References")),
        },
        "edits": [],
    }


def merge_provenance(rows):
    """Preserve every original evidence record when parallel rows consolidate."""
    values = [provenance(row) for row in rows]
    result = copy.deepcopy(values[0])
    originals, edits = [], []
    for value in values:
        for original in value.get("originals", [value["original"]]):
            if original not in originals:
                originals.append(original)
        for edit in value["edits"]:
            if edit not in edits:
                edits.append(edit)
    if len(originals) > 1:
        result["originals"] = originals
    result["edits"] = edits
    return result


def delta(before, after, keys):
    """Exact records; ambiguous parallel identities are additions/removals."""
    old, new = {}, {}
    for target, frame in ((old, before), (new, after)):
        for row in records(frame):
            target.setdefault(tuple(row.get(k) for k in keys), []).append(row)
    added, removed, changed = [], [], []
    for key in sorted(old.keys() | new.keys(), key=str):
        left, right = old.get(key, []), new.get(key, [])
        left_only = [r for r in left if r not in right]
        right_only = [r for r in right if r not in left]
        if len(left) == len(right) == 1 and left_only and right_only:
            changed.append({"before": left[0], "after": right[0]})
        else:
            removed.extend(left_only)
            added.extend(right_only)
    return {"added": added, "removed": removed, "changed": changed}


def bound_parameters(method, *args, **kwargs):
    bound = inspect.signature(method).bind(*args, **kwargs)
    bound.apply_defaults()
    return {
        k: json_value(v)
        for k, v in bound.arguments.items()
        if k not in {"self", "network"} and not k.startswith("_")
    }


def construction_operation(method):
    """Record nested construction calls and their actual effective inputs."""

    @wraps(method)
    def wrapped(network, *args, **kwargs):
        context = getattr(network, "_transition_context", None)
        if context is None:
            return method(network, *args, **kwargs)
        operation = {
            "method": method.__name__,
            "parameters": bound_parameters(method, network, *args, **kwargs),
        }
        context["operations"].append(operation)
        context["stack"].append(operation)
        try:
            return method(network, *args, **kwargs)
        finally:
            context["stack"].pop()

    return wrapped


def effective_parameters(network, **values):
    context = getattr(network, "_transition_context", None)
    if context and context["stack"]:
        context["stack"][-1]["parameters"].update(json_value(values))


def inclusion(network, source, target, reason):
    context = getattr(network, "_transition_context", None)
    if context is None:
        return
    item = {
        "source": source,
        "target": target,
        **json_value(reason),
        "operation_indices": [
            next(
                i
                for i, candidate in enumerate(context["operations"])
                if candidate is operation
            )
            for operation in context["stack"]
        ],
    }
    # Keep multiple actual paths, but avoid identical repeated records.
    if item not in context["reasons"]:
        context["reasons"].append(item)


def history_payload(network):
    return {
        "format": "neko-history",
        "schema_version": 1,
        "database_name": network.database_name,
        "current_state_id": network.current_state_id,
        "root_state_id": network.root_state_id,
        "next_state_id": network._state_counter,
        "max_history": network._max_history,
        "pruned_state_ids": list(network._pruned_state_ids),
        "history_enabled": network._history_enabled,
        "working_state_matches_current": working_state_matches_current(network),
        "working_state": {
            "nodes": frame_payload(network.nodes),
            "edges": frame_payload(network.edges),
            "initial_nodes": json_value(network.initial_nodes),
        },
        "states": [
            {
                "id": sid,
                "parent_ids": list(state.parent_ids),
                "children_ids": list(state.children_ids),
                "nodes": frame_payload(state.nodes),
                "edges": frame_payload(state.edges),
                "metadata": json_value(state.metadata),
            }
            for sid in network._state_log
            for state in [network._states[sid]]
        ],
    }


def working_state_matches_current(network):
    state = network._states.get(network.current_state_id)
    return state is not None and (
        frame_fingerprint(state.nodes) == frame_fingerprint(network.nodes)
        and frame_fingerprint(state.edges) == frame_fingerprint(network.edges)
        and state.metadata.get("initial_nodes")
        == json_value(network.initial_nodes)
    )


def validate_snapshot(nodes, edges):
    """Reject broken snapshots before publishing a recovered Network."""
    if not {"Uniprot", "Genesymbol", "Type"} <= set(nodes.columns):
        raise ValueError("History nodes are missing required columns.")
    if not {"source", "target", "Type", "Effect", "References"} <= set(
        edges.columns
    ):
        raise ValueError("History edges are missing required columns.")
    identifiers = set(nodes.Uniprot.dropna())
    if any(not isinstance(value, str) or not value for value in identifiers):
        raise ValueError("History contains invalid node identifiers.")
    for row in records(edges):
        if row["source"] not in identifiers or row["target"] not in identifiers:
            raise ValueError(
                "History contains an edge with a missing endpoint."
            )
        value = row.get("Provenance")
        if value is not None:
            validate_provenance(value)


def validate_provenance(value):
    if (
        not isinstance(value, dict)
        or not isinstance(value.get("origin"), str)
        or value.get("origin")
        not in {"manual", "resource", "import", "unknown"}
        or not isinstance(value.get("original"), dict)
        or not isinstance(value.get("edits"), list)
        or any(not isinstance(edit, dict) for edit in value["edits"])
    ):
        raise ValueError("Malformed NeKo interaction provenance.")
    json_value(value)


def write_json(path, payload):
    path = Path(path)
    text = json.dumps(json_value(payload), indent=2, allow_nan=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text + "\n", encoding="utf-8")
    return str(path)


def load_history(cls, path, resources=None):
    from .network_state import NetworkState

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if (
        isinstance(payload, dict)
        and payload.get("format") == "neko-export"
        and payload.get("schema_version") == 1
    ):
        payload = payload["history"]
    if (
        not isinstance(payload, dict)
        or payload.get("format") != "neko-history"
        or payload.get("schema_version") != 1
    ):
        raise ValueError("Unsupported NeKo history format or schema version.")
    states = {}
    for item in payload["states"]:
        sid = item["id"]
        if type(sid) is not int or sid < 0 or sid in states:
            raise ValueError("Invalid or duplicate history state id.")
        state = NetworkState(
            frame_from_payload(item["nodes"]),
            frame_from_payload(item["edges"]),
            item["metadata"],
            sid,
            item["parent_ids"],
        )
        validate_snapshot(state.nodes, state.edges)
        if state.metadata.get("database_name") != payload["database_name"]:
            raise ValueError(
                "History states have inconsistent database identities."
            )
        state.children_ids = item["children_ids"]
        states[sid] = state
    if (
        not states
        or payload["current_state_id"] not in states
        or payload["root_state_id"] not in states
    ):
        raise ValueError("History must contain its root and current state.")
    import networkx as nx

    graph = nx.DiGraph()
    graph.add_nodes_from(states)
    for sid, state in states.items():
        for parent in state.parent_ids:
            if parent not in states or sid not in states[parent].children_ids:
                raise ValueError(
                    "Inconsistent history parent/child relationship."
                )
            graph.add_edge(parent, sid)
        for child in state.children_ids:
            if child not in states or sid not in states[child].parent_ids:
                raise ValueError(
                    "Inconsistent history child/parent relationship."
                )
    if not nx.is_directed_acyclic_graph(graph):
        raise ValueError("History contains a cycle.")
    root = payload["root_state_id"]
    if states[root].parent_ids or nx.descendants(graph, root) | {root} != set(
        states
    ):
        raise ValueError("History contains states disconnected from its root.")
    if type(payload["next_state_id"]) is not int or payload[
        "next_state_id"
    ] <= max(states):
        raise ValueError("Invalid next history state id.")
    # An empty resource enables offline restoration without acquiring a DB.
    empty = pd.DataFrame(
        columns=[
            "source",
            "target",
            "is_directed",
            "is_stimulation",
            "is_inhibition",
            "form_complex",
        ]
    )
    working = payload["working_state"]
    working_nodes = frame_from_payload(working["nodes"])
    working_edges = frame_from_payload(working["edges"])
    validate_snapshot(working_nodes, working_edges)
    net = cls(
        resources=resources if resources is not None else empty,
        database_name=payload["database_name"],
    )
    net._resources_available = resources is not None
    net._states = states
    net._state_log = list(states)
    net._state_metadata = {sid: state.metadata for sid, state in states.items()}
    net._root_state_id = payload["root_state_id"]
    net._state_counter = payload["next_state_id"]
    net._max_history = payload["max_history"]
    net._pruned_state_ids = payload["pruned_state_ids"]
    net._history_enabled = payload["history_enabled"]
    net.checkout(payload["current_state_id"])
    net.nodes = working_nodes
    net.edges = working_edges
    net.initial_nodes = working["initial_nodes"]
    net.sync_nodes_from_df()
    net.sync_edges_from_df()
    return net
