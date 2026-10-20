"""Bounded history work and unchanged scientific/provenance semantics."""

import copy
import json
from unittest.mock import patch

import pandas as pd
import pytest

from neko.core.network import Network
from neko.core import history
from neko._outputs.exports import Exports


@pytest.fixture(autouse=True)
def offline_identifiers(monkeypatch):
    monkeypatch.setattr(
        "neko.core.network.mapping_node_identifier",
        lambda value: [None, value, value],
    )
    monkeypatch.setattr(
        "neko.core.network.check_gene_list_format",
        lambda values: False,
    )


def ring(size, seeds=None):
    resources = pd.DataFrame(
        [
            {
                "source": f"N{i}",
                "target": f"N{(i + 1) % size}",
                "is_directed": True,
                "is_stimulation": True,
                "is_inhibition": False,
                "form_complex": False,
                "references": "123456",
            }
            for i in range(size)
        ]
    )
    return Network(
        initial_nodes=seeds or [f"N{i}" for i in range(size)],
        resources=resources,
    )


def forbidden(*args, **kwargs):
    raise AssertionError("Unexpected full-history traversal or copy")


class Uncopyable:
    def __deepcopy__(self, memo):
        forbidden()


def test_status_uses_only_current_snapshot_and_detects_nested_edits(
    monkeypatch,
):
    network = ring(5)
    network.connect_nodes()
    network.modify_node_name("N0", "Renamed")
    assert network.state_count == 3
    network.set_max_history(2)
    assert network.state_count == 2
    monkeypatch.setattr(network, "history_dict", forbidden)
    monkeypatch.setattr(history, "history_payload", forbidden)
    # An old snapshot cannot be consulted to answer current working status.
    network._states[network.root_state_id].edges = Uncopyable()
    assert network.history_status() == {
        "pruned_state_count": 1,
        "working_state_matches_current": True,
    }
    network.edges.at[0, "Provenance"]["original"]["references"] = "654321"
    assert network.history_status()["working_state_matches_current"] is False
    assert (
        network.get_state()["edges"]["records"][0]["Provenance"]["original"][
            "references"
        ]
        == "123456"
    )


def test_pagination_copies_only_selected_states_and_summary_excludes_evidence():
    network = ring(5)
    network.connect_nodes()
    network.modify_node_name("N0", "Renamed")
    full = network.list_states()
    assert network.list_states(offset=1, limit=1) == full[1:2]
    assert network.list_states(offset=999) == []
    assert network.list_states(limit=0) == []
    network._state_metadata[0]["unselected"] = Uncopyable()
    assert network.list_states(offset=1, limit=1) == full[1:2]
    network._state_metadata[1]["transition"]["operations"] = Uncopyable()
    summary = network.list_states(offset=1, limit=1, detail="summary")[0]
    assert summary["metadata"]["transition"]["counts"]["edges"]["added"] == 5
    assert "operations" not in summary["metadata"]["transition"]
    assert "pruned_state_ids" not in summary
    summary["parent_ids"].clear()
    summary["metadata"]["transition"]["counts"]["edges"]["added"] = 999
    assert network.list_states(offset=1, limit=1, detail="summary")[0][
        "parent_ids"
    ] == [0]
    assert len(network._states[1].metadata["transition"]["edges"]["added"]) == 5


@pytest.mark.parametrize(
    "kwargs",
    [
        {"offset": -1},
        {"offset": True},
        {"offset": 1.5},
        {"limit": -1},
        {"limit": False},
        {"limit": "5"},
        {"detail": "brief"},
    ],
)
def test_invalid_page_arguments(kwargs):
    with pytest.raises(ValueError):
        ring(3).list_states(**kwargs)


def test_summary_does_not_invent_change_counts_for_legacy_snapshots():
    network = ring(3)
    network._state_metadata[0]["transition"] = {
        "explanation_status": "unavailable",
    }
    transition = network.list_states(detail="summary")[0]["metadata"]["transition"]
    assert transition["counts"] == {"nodes": None, "edges": None}


def test_mutation_marks_preexisting_nested_edit_and_keeps_past_evidence():
    network = ring(5)
    network.connect_nodes()
    old = network.get_state()
    network.edges.at[0, "Provenance"]["original"]["references"] = "654321"
    network.modify_node_name("N0", "Renamed")
    assert network.get_state(old["id"])["edges"] == old["edges"]
    assert (
        network.get_state()["metadata"]["transition"]["explanation_status"]
        == "partial"
    )
    assert network.history_status()["working_state_matches_current"] is True


def test_export_capture_is_independent_and_never_copies_network(
    tmp_path, monkeypatch
):
    network = ring(5)
    network.connect_nodes()
    captured = network.history_dict()
    monkeypatch.setattr(network, "copy", forbidden)
    exporter = Exports(network)
    network.edges.at[0, "Provenance"]["original"]["references"] = "654321"
    network.modify_node_name("N0", "Renamed")
    network.set_max_history(2)
    monkeypatch.setattr(history, "history_payload", forbidden)
    for suffix in ("first", "second"):
        path = tmp_path / f"{suffix}.sif"
        exporter.export_sif(path)
        payload = json.loads((tmp_path / f"{suffix}.sif.neko.json").read_text())
        assert payload["history"] == captured
        assert (
            payload["export"]["interactions"][0]["Provenance"]["original"][
                "references"
            ]
            == "123456"
        )
        assert Network.from_history(str(path) + ".neko.json").history_status()[
            "working_state_matches_current"
        ]


def test_compact_and_pretty_history_round_trip_and_failed_write_is_atomic(
    tmp_path,
):
    network = ring(5)
    network.connect_nodes()
    compact, pretty = tmp_path / "compact.json", tmp_path / "pretty.json"
    network.export_history(compact)
    network.export_history(pretty, pretty=True)
    assert compact.stat().st_size < pretty.stat().st_size
    assert json.loads(compact.read_text()) == json.loads(pretty.read_text())
    assert (
        Network.from_history(compact).history_dict() == network.history_dict()
    )
    before = compact.read_bytes()
    with pytest.raises(ValueError):
        history.write_json(compact, {"invalid": float("inf")})
    assert compact.read_bytes() == before
    assert set(tmp_path.iterdir()) == {compact, pretty}


@pytest.mark.parametrize("size", [8, 40])
@pytest.mark.parametrize(
    "strategy",
    [
        "connect_nodes",
        "complete",
        "radial",
        "upstream",
        "subgroup",
        "atopo",
    ],
)
def test_strategy_captures_once_without_per_edge_fingerprints(size, strategy):
    seeds = None if strategy == "connect_nodes" else ["N0", f"N{size - 1}"]
    network = ring(size, seeds)
    untracked = network.copy()
    untracked.set_history_tracking(False)

    def grow(net):
        if strategy == "connect_nodes":
            net.connect_nodes()
        elif strategy == "complete":
            net.complete_connection(
                maxlen=size,
                path_policy="all_shortest",
                reuse_policy="discovered_paths",
            )
        elif strategy == "radial":
            net.connect_network_radially(max_len=size, direction="OUT")
        elif strategy == "upstream":
            net.connect_to_upstream_nodes(depth=size, rank=2)
        elif strategy == "subgroup":
            net.connect_subgroup(seeds, maxlen=size)
        else:
            net.connect_as_atopo(strategy="complete", max_len=size)

    with patch.object(
        history, "frame_fingerprint", wraps=history.frame_fingerprint
    ) as fingerprint:
        with patch.object(
            network, "save_state", wraps=network.save_state
        ) as save:
            with patch.object(
                history, "history_payload", side_effect=forbidden
            ):
                grow(network)
            assert save.call_count == 1
        # Two working tables before/after, and the current snapshot's two tables.
        assert fingerprint.call_count == 6
    grow(untracked)
    assert history.records(network.edges) == history.records(untracked.edges)
    transition = network.get_state()["metadata"]["transition"]
    assert transition["edges"]["added"]
    assert transition["operations"]
    for reason in transition["inclusion_reasons"]:
        assert reason["kind"] != "unavailable"
        assert all(
            0 <= i < len(transition["operations"])
            for i in reason["operation_indices"]
        )
    assert network._transition_context is None


def test_unique_pairs_never_serialize_provenance_for_deduplication(monkeypatch):
    network = ring(5)
    network.connect_nodes()
    expected = network.edges.reset_index(drop=True)
    monkeypatch.setattr(history, "records", forbidden)
    pd.testing.assert_frame_equal(history.unique_edges(network.edges), expected)


def test_duplicate_pairs_preserve_distinct_evidence_effects_and_order():
    network = ring(5)
    network.connect_nodes()
    rows = history.records(network.edges)
    different_evidence = copy.deepcopy(rows[0])
    different_evidence["Provenance"]["original"]["references"] = "654321"
    different_effect = copy.deepcopy(rows[0])
    different_effect["Effect"] = "inhibition"
    table = pd.DataFrame(
        [
            rows[0],
            rows[1],
            rows[0],
            different_evidence,
            different_effect,
            rows[2],
        ],
        index=[5] * 6,
    )
    assert history.records(history.unique_edges(table)) == [
        rows[0],
        rows[1],
        different_evidence,
        different_effect,
        rows[2],
    ]


def test_reason_deduplication_does_not_scan_previous_reasons():
    class NoContains(list):
        def __contains__(self, value):
            raise AssertionError("Linear explanation scan")

    network = ring(3)
    operation = {"method": "grow", "parameters": {}}
    network._transition_context = {
        "operations": [operation],
        "stack": [operation],
        "reasons": NoContains(),
    }
    for _ in range(2):
        for middle in ("X", "Y"):
            history.inclusion(
                network,
                "N0",
                "N1",
                {
                    "kind": "selected_path",
                    "path": ["N0", middle, "N1"],
                },
            )
    assert len(network._transition_context["reasons"]) == 2
    assert all(
        r["operation_indices"] == [0]
        for r in network._transition_context["reasons"]
    )


@pytest.mark.parametrize(
    "effect, expected",
    [
        ("activate", "stimulation"),
        ("inhibit", "inhibition"),
        ("both", "bimodal"),
        ("form_complex", "form complex"),
        (None, "undefined"),
        (" custom ", "custom"),
    ],
)
def test_single_edge_consolidation_preserves_normalization_and_evidence(
    effect, expected
):
    from neko.core.tools import consolidate_edges

    frame = pd.DataFrame(
        [
            {
                "source": "N0",
                "target": "N1",
                "Effect": effect,
                "Type": " A; A; B ",
                "References": "123; 456; 123",
            }
        ]
    )
    frame["Provenance"] = [history.provenance(frame.iloc[0])]
    result = consolidate_edges(frame)
    assert result.iloc[0]["Effect"] == expected
    assert result.iloc[0]["Type"] == "A; B"
    assert result.iloc[0]["References"] == "123; 456"
    assert result.iloc[0]["Provenance"] == frame.iloc[0]["Provenance"]
    result.at[0, "Provenance"]["original"]["references"] = "changed"
    assert (
        frame.at[0, "Provenance"]["original"]["references"] == "123; 456; 123"
    )
