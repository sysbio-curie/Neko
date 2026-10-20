"""Validated, atomic manual interaction curation owned by NeKo."""

import copy
import re
from typing import Any, Optional
from typing_extensions import Literal, NotRequired, TypedDict

import pandas as pd

from .history import (
    copy_frame,
    json_value,
    provenance,
    records,
    working_state_matches_current,
)

EFFECTS = {"stimulation", "inhibition", "bimodal", "form complex", "undefined"}


class InteractionEvidence(TypedDict):
    summary: str
    pmids: NotRequired[list[str]]
    citations: NotRequired[list[str]]


class InteractionEdit(TypedDict):
    operation: Literal["add", "set_effect"]
    source: str
    target: str
    effect: Literal[
        "stimulation", "inhibition", "bimodal", "form complex", "undefined"
    ]
    rationale: str
    evidence: InteractionEvidence


class InteractionChange(TypedDict):
    operation: Literal["add", "set_effect"]
    source: str
    target: str
    before: Optional[dict[str, Any]]
    after: dict[str, Any]


class InteractionEditResult(TypedDict):
    changed: bool
    previous_state_id: int
    state_id: int
    changes: list[InteractionChange]


class InteractionEditError(ValueError):
    """Machine-readable failure; no part of a rejected batch is applied."""

    def __init__(self, code, message, edit_index=None, details=None):
        super().__init__(message)
        self.code = code
        self.edit_index = edit_index
        self.details = details or {}


def _fail(code, message, index=None, **details):
    raise InteractionEditError(code, message, index, details)


def _endpoint(network, identifier, index):
    if not isinstance(identifier, str) or not identifier.strip():
        _fail("invalid_request", "Endpoints must be nonempty strings.", index)
    matches = network.nodes[
        (network.nodes["Uniprot"] == identifier)
        | (network.nodes["Genesymbol"] == identifier)
    ]
    ids = {
        row["Uniprot"] if pd.notna(row["Uniprot"]) else row["Genesymbol"]
        for _, row in matches.iterrows()
    }
    if not ids:
        _fail(
            "endpoint_not_found",
            "Endpoint is not in this network.",
            index,
            identifier=identifier,
        )
    if len(ids) != 1:
        _fail(
            "ambiguous_endpoint",
            "Endpoint identifies multiple nodes.",
            index,
            identifier=identifier,
        )
    return next(iter(ids))


def _evidence(value, index):
    if not isinstance(value, dict) or set(value) - {
        "summary",
        "pmids",
        "citations",
    }:
        _fail(
            "invalid_evidence",
            "Evidence must contain summary, pmids and/or citations.",
            index,
        )
    summary = value.get("summary")
    if not isinstance(summary, str) or not summary.strip():
        _fail(
            "invalid_evidence",
            "A nonempty evidence summary is required.",
            index,
        )
    pmids, citations = value.get("pmids", []), value.get("citations", [])
    if not isinstance(pmids, list) or any(
        not isinstance(p, str) or not re.fullmatch(r"[1-9][0-9]*", p)
        for p in pmids
    ):
        _fail(
            "invalid_evidence",
            "PMIDs must be positive decimal strings, without a prefix.",
            index,
        )
    if not isinstance(citations, list) or any(
        not isinstance(c, str) or not c.strip() for c in citations
    ):
        _fail("invalid_evidence", "Citations must be nonempty strings.", index)
    if not pmids and not citations:
        _fail(
            "invalid_evidence", "Supply a PMID or a supporting citation.", index
        )
    return {
        "summary": summary.strip(),
        "pmids": list(dict.fromkeys(pmids)),
        "citations": list(dict.fromkeys(c.strip() for c in citations)),
    }


def edit_interactions(
    network,
    edits: list[InteractionEdit],
    expected_state_id: Optional[int] = None,
) -> InteractionEditResult:
    if not network._history_enabled:
        _fail("history_disabled", "Manual curation requires history tracking.")
    if expected_state_id is not None:
        if type(expected_state_id) is not int:
            _fail("invalid_request", "expected_state_id must be an integer.")
        if expected_state_id != network.current_state_id:
            _fail(
                "stale_state",
                "The current state differs from the expected state.",
                expected=expected_state_id,
                actual=network.current_state_id,
            )
    if not isinstance(edits, list) or not edits:
        _fail("invalid_request", "edits must be a nonempty list.")
    frame = copy_frame(network.edges)
    if "Provenance" not in frame:
        frame["Provenance"] = [provenance(row) for row in records(frame)]
    frame = frame.reset_index(drop=True)
    changes, seen = [], set()
    next_id = network._state_counter
    for index, edit in enumerate(edits):
        required = {
            "operation",
            "source",
            "target",
            "effect",
            "rationale",
            "evidence",
        }
        if not isinstance(edit, dict) or set(edit) != required:
            _fail(
                "invalid_request",
                "Each edit requires exactly operation, source, target, effect, rationale and evidence.",
                index,
            )
        operation = edit["operation"]
        if operation not in ("add", "set_effect"):
            _fail(
                "invalid_request", "operation must be add or set_effect.", index
            )
        if not isinstance(edit["effect"], str) or edit["effect"] not in EFFECTS:
            _fail(
                "invalid_effect",
                "Unsupported interaction effect.",
                index,
                allowed=sorted(EFFECTS),
            )
        rationale = edit["rationale"]
        if not isinstance(rationale, str) or not rationale.strip():
            _fail("invalid_request", "A nonempty rationale is required.", index)
        evidence = _evidence(edit["evidence"], index)
        source, target = (
            _endpoint(network, edit[k], index) for k in ("source", "target")
        )
        if (source, target) in seen:
            _fail(
                "invalid_request",
                "A batch may edit each directed pair only once.",
                index,
            )
        seen.add((source, target))
        matches = frame[(frame.source == source) & (frame.target == target)]
        if operation == "add" and not matches.empty:
            _fail(
                "interaction_exists",
                "Use set_effect to edit an existing interaction.",
                index,
            )
        if operation == "set_effect" and matches.empty:
            _fail(
                "interaction_not_found",
                "No interaction exists for this directed pair.",
                index,
            )
        if len(matches) > 1:
            _fail(
                "ambiguous_interaction",
                "Multiple interactions exist for this directed pair.",
                index,
            )
        before = records(matches)[0] if not matches.empty else None
        if before and before["Effect"] == edit["effect"]:
            continue
        after = (
            copy.deepcopy(before)
            if before
            else {
                "source": source,
                "target": target,
                "Type": None,
                "Effect": edit["effect"],
                "References": ";".join(evidence["pmids"]) or None,
            }
        )
        after["Effect"] = edit["effect"]
        prov = provenance(
            before if before else after,
            origin="manual" if before is None else "unknown",
        )
        prov["edits"].append(
            {
                "operation": operation,
                "origin": "manual",
                "state_id": next_id,
                "before_effect": before["Effect"] if before else None,
                "after_effect": edit["effect"],
                "rationale": rationale.strip(),
                "evidence": evidence,
            }
        )
        after["Provenance"] = prov
        if before:
            for key, value in after.items():
                frame.at[matches.index[0], key] = value
        else:
            frame = pd.concat([frame, pd.DataFrame([after])], ignore_index=True)
        changes.append(
            {
                "operation": operation,
                "source": source,
                "target": target,
                "before": before,
                "after": after,
            }
        )
    previous = network.current_state_id
    if changes:
        # Only publish after the complete batch is validated and JSON-safe.
        changes = json_value(changes)
        # Stage only mutable working/history bookkeeping. Resource tables and
        # connection indexes can be large and are not touched by curation.
        candidate = copy.copy(network)
        candidate._states = {
            sid: copy.copy(state) for sid, state in network._states.items()
        }
        for state in candidate._states.values():
            state.parent_ids = list(state.parent_ids)
            state.children_ids = list(state.children_ids)
            state.metadata = copy.deepcopy(state.metadata)
        candidate._state_metadata = {
            sid: state.metadata for sid, state in candidate._states.items()
        }
        candidate._state_log = list(network._state_log)
        candidate._pruned_state_ids = list(network._pruned_state_ids)
        candidate.edges = frame
        candidate.sync_nodes_from_df()
        candidate.sync_edges_from_df()
        candidate.save_state(
            metadata={
                "method": "edit_interactions",
                "label": "Manual interaction curation",
                "manual_changes": changes,
                "untracked_changes": not working_state_matches_current(network),
            }
        )
        for name in (
            "edges",
            "_node_objs",
            "_edge_objs",
            "_states",
            "_state_metadata",
            "_state_log",
            "_state_counter",
            "_current_state_id",
            "_root_state_id",
            "_pruned_state_ids",
        ):
            setattr(network, name, getattr(candidate, name))
    return {
        "changed": bool(changes),
        "previous_state_id": previous,
        "state_id": network.current_state_id,
        "changes": copy.deepcopy(changes),
    }
