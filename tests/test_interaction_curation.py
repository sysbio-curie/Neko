"""Offline scientific contract for W2 interaction curation and state recovery."""

import copy
import json
import pickle

import numpy as np
import pandas as pd
import pytest

from neko.core.network import Network, InteractionEditError
from neko.core.history import history_payload, records
from neko._outputs.exports import Exports
from neko.inputs import _universe


@pytest.fixture(autouse=True)
def identifiers(monkeypatch):
    import neko.core.network as module

    monkeypatch.setattr(
        module,
        "mapping_node_identifier",
        lambda value: [
            None,
            value.removeprefix("ID:"),
            value if value.startswith("ID:") else "ID:" + value,
        ],
    )
    monkeypatch.setattr(module, "check_gene_list_format", lambda values: True)


@pytest.fixture
def resources():
    return pd.DataFrame(
        [
            {
                "source": "ID:" + a,
                "target": "ID:" + b,
                "type": "regulation",
                "is_stimulation": True,
                "is_inhibition": False,
                "is_directed": True,
                "form_complex": False,
                "references": str(index + 100),
            }
            for index, (a, b) in enumerate(
                [("A", "B"), ("B", "C"), ("A", "D"), ("D", "C"), ("C", "E")]
            )
        ]
    )


@pytest.fixture
def network(resources):
    return Network(
        initial_nodes=["ID:A", "ID:B", "ID:C"],
        resources=resources,
        database_name="SIGNOR",
    )


def edit(
    operation="add", source="A", target="C", effect="inhibition", pmids=None
):
    return {
        "operation": operation,
        "source": source,
        "target": target,
        "effect": effect,
        "rationale": "Explicit literature-supported model correction.",
        "evidence": {
            "summary": "The cited experiment supports this effect.",
            "pmids": ["123456"] if pmids is None else pmids,
            "citations": ["Example publication, section 2"],
        },
    }


def assert_objects(network):
    assert {n.id for n in network.nodes_as_objects()} == set(
        network.nodes.Uniprot
    )
    expected = {
        (r["source"], r["target"], r["Type"], r["Effect"]): r
        for r in records(network.edges)
    }
    actual = {
        (e.source, e.target, e.interaction_type, e.metadata["Effect"]): e
        for e in network.edges_as_objects()
    }
    assert actual.keys() == expected.keys()
    for key, obj in actual.items():
        row = expected[key]
        assert obj.evidence == row["References"]
        assert obj.metadata["Provenance"] == row["Provenance"]


def test_manual_absent_database_edge_requires_no_lookup_or_growth(
    network, monkeypatch
):
    import neko.core.network as module

    before_nodes = records(network.nodes)
    monkeypatch.setattr(
        module,
        "mapping_node_identifier",
        lambda *_: pytest.fail("external lookup"),
    )
    monkeypatch.setattr(
        network, "add_node", lambda *_: pytest.fail("implicit node creation")
    )
    result = network.edit_interactions(
        [edit()], expected_state_id=network.current_state_id
    )
    assert (
        result["changed"] and result["state_id"] != result["previous_state_id"]
    )
    assert records(network.nodes) == before_nodes
    edge = result["changes"][0]["after"]
    assert (edge["source"], edge["target"]) == ("ID:A", "ID:C")
    assert edge["Provenance"]["origin"] == "manual"
    assert edge["Provenance"]["original"]["database_name"] is None
    assert edge["Provenance"]["edits"][0]["evidence"]["pmids"] == ["123456"]
    assert edge["Provenance"]["edits"][0]["state_id"] == result["state_id"]
    assert edge["References"] == "123456"
    assert_objects(network)


@pytest.mark.parametrize(
    "effect",
    ["stimulation", "inhibition", "bimodal", "form complex", "undefined"],
)
def test_supported_effects(network, effect):
    assert (
        network.edit_interactions([edit(effect=effect)])["changes"][0]["after"][
            "Effect"
        ]
        == effect
    )


def test_no_pmid_does_not_invent_database_evidence(network):
    result = network.edit_interactions([edit(pmids=[])])
    edge = result["changes"][0]["after"]
    assert edge["References"] is None
    assert edge["Provenance"]["edits"][0]["evidence"]["pmids"] == []
    assert edge["Provenance"]["origin"] == "manual"


@pytest.mark.parametrize(
    ("replacement", "code"),
    [
        ({"target": "MISSING"}, "endpoint_not_found"),
        ({"source": ""}, "invalid_request"),
        ({"operation": "remove"}, "invalid_request"),
        ({"operation": "set_effect"}, "interaction_not_found"),
        ({"effect": "+"}, "invalid_effect"),
        ({"effect": []}, "invalid_effect"),
        ({"rationale": " "}, "invalid_request"),
        (
            {
                "evidence": {
                    "summary": "unsupported",
                    "pmids": [],
                    "citations": [],
                }
            },
            "invalid_evidence",
        ),
        (
            {"evidence": {"summary": "test", "pmids": ["PMID:123"]}},
            "invalid_evidence",
        ),
        ({"evidence": {"summary": "test", "pmids": [123]}}, "invalid_evidence"),
        ({"evidence": {"summary": "", "pmids": ["123"]}}, "invalid_evidence"),
    ],
)
def test_invalid_batch_is_atomic(network, replacement, code):
    before = history_payload(network)
    invalid = edit(source="B", target="A")
    invalid.update(replacement)
    with pytest.raises(InteractionEditError) as exc:
        network.edit_interactions([edit(), invalid])
    assert exc.value.code == code
    assert exc.value.edit_index == 1
    assert history_payload(network) == before
    assert network.edges_as_objects() == set()


def test_duplicate_add_and_duplicate_pair_fail_without_overwrite(network):
    with pytest.raises(InteractionEditError, match="only once"):
        network.edit_interactions([edit(), edit()])
    network.edit_interactions([edit()])
    before = history_payload(network)
    with pytest.raises(InteractionEditError) as exc:
        network.edit_interactions([edit(effect="stimulation")])
    assert exc.value.code == "interaction_exists"
    assert history_payload(network) == before


def test_ambiguous_endpoint_and_interaction(network):
    network.nodes.loc[network.nodes.Uniprot == "ID:B", "Genesymbol"] = "A"
    with pytest.raises(InteractionEditError) as exc:
        network.edit_interactions([edit()])
    assert exc.value.code == "ambiguous_endpoint"
    network.nodes.loc[network.nodes.Uniprot == "ID:B", "Genesymbol"] = "B"
    network.connect_nodes()
    extra = network.edges.iloc[[0]].copy()
    extra["Effect"] = "form complex"
    network.edges = pd.concat([network.edges, extra], ignore_index=True)
    with pytest.raises(InteractionEditError) as exc:
        network.edit_interactions([edit("set_effect", source="A", target="B")])
    assert exc.value.code == "ambiguous_interaction"


def test_stale_disabled_and_noop(network):
    initial = network.current_state_id
    network.edit_interactions([edit()])
    with pytest.raises(InteractionEditError) as exc:
        network.edit_interactions(
            [edit("set_effect")], expected_state_id=initial
        )
    assert exc.value.code == "stale_state"
    with pytest.raises(InteractionEditError) as exc:
        network.edit_interactions([edit("set_effect")], expected_state_id=True)
    assert exc.value.code == "invalid_request"
    before = history_payload(network)
    result = network.edit_interactions([edit("set_effect")])
    assert result == {
        "changed": False,
        "changes": [],
        "state_id": network.current_state_id,
        "previous_state_id": network.current_state_id,
    }
    assert history_payload(network) == before
    with network.suspend_history(), pytest.raises(InteractionEditError) as exc:
        network.edit_interactions([edit("set_effect", effect="stimulation")])
    assert exc.value.code == "history_disabled"


def test_effect_change_retains_database_evidence_and_survives_growth(network):
    network.connect_nodes()
    before = records(network.edges)[0]
    result = network.edit_interactions([edit("set_effect", target="B")])
    after = result["changes"][0]["after"]
    assert after["References"] == before["References"] == "100"
    assert after["Provenance"]["origin"] == "resource"
    assert after["Provenance"]["original"]["effect"] == "stimulation"
    assert after["Provenance"]["original"]["database_name"] == "SIGNOR"
    assert after["Provenance"]["edits"][0]["origin"] == "manual"
    assert after["Provenance"]["edits"][0]["evidence"]["pmids"] == ["123456"]
    network.connect_nodes()
    row = network.edges.loc[network.edges.source == "ID:A"].iloc[0]
    assert row.Effect == "inhibition"
    assert row.Provenance["edits"] == after["Provenance"]["edits"]
    assert_objects(network)


def test_all_errors_including_snapshot_failure_are_atomic(network, monkeypatch):
    before = history_payload(network)
    monkeypatch.setattr(
        Network,
        "save_state",
        lambda *_args, **_kw: (_ for _ in ()).throw(
            RuntimeError("disk-like failure")
        ),
    )
    with pytest.raises(RuntimeError):
        network.edit_interactions([edit()])
    assert history_payload(network) == before


def test_legacy_rejected_endpoint_has_no_partial_node_or_dangling_edge(network):
    before = history_payload(network)
    row = pd.DataFrame(
        [
            {
                "source": "ID:D",
                "target": "ID:MISSING",
                "is_stimulation": True,
                "is_inhibition": False,
            }
        ]
    )
    with pytest.raises(ValueError, match="rejected"):
        network.add_edge(row)
    assert history_payload(network) == before
    assert set(n.id for n in network.nodes_as_objects()) == set(
        network.nodes.Uniprot
    )


def test_branches_metadata_isolation_and_comparison(network):
    network.connect_nodes()
    parent = network.current_state_id
    first = network.edit_interactions([edit()])["state_id"]
    network.checkout(parent)
    assert_objects(network)
    second = network.edit_interactions([edit(effect="stimulation")])["state_id"]
    difference = network.compare_states(first, second)
    changed = difference["edge_changes"]["changed"]
    assert len(changed) == 1
    assert changed[0]["before"]["Effect"] == "inhibition"
    assert changed[0]["after"]["Effect"] == "stimulation"
    assert set(network.history_graph().successors(parent)) == {first, second}
    network.checkout(first)
    assert_objects(network)
    network.edges.iloc[-1].Provenance["edits"][0][
        "rationale"
    ] = "mutated live object"
    assert (
        network._states[first]
        .edges.iloc[-1]
        .Provenance["edits"][0]["rationale"]
        != "mutated live object"
    )
    returned = network.list_states()
    returned[-1]["metadata"]["manual_changes"].clear()
    assert network._states[second].metadata["manual_changes"]
    network.checkout(second)
    clone = network.copy()
    clone.edges.iloc[-1].Provenance["edits"].clear()
    assert network.edges.iloc[-1].Provenance["edits"]
    network.sync_edges_to_df()
    assert_objects(network)


@pytest.mark.parametrize(
    "policy", ["one_shortest", "all_shortest", "all_bounded"]
)
def test_actual_path_attribution_and_effective_defaults(resources, policy):
    network = Network(initial_nodes=["ID:A", "ID:C"], resources=resources)
    network.complete_connection(
        path_policy=policy, reuse_policy="none", maxlen=2
    )
    transition = network.list_states()[-1]["metadata"]["transition"]
    operation = transition["operations"][0]
    assert operation["parameters"]["path_policy"] == policy
    assert operation["parameters"]["only_signed"] is False
    assert operation["parameters"]["consensus"] is False
    assert "algorithm" not in operation["parameters"]
    reasons = transition["inclusion_reasons"]
    for edge in transition["edges"]["added"]:
        assert any(
            (r["source"], r["target"]) == (edge["source"], edge["target"])
            and r["kind"] != "unavailable"
            for r in reasons
        )
    if policy == "all_shortest":
        union = [r for r in reasons if r["kind"] == "shortest_path_union"]
        assert union and all(r["seed_pair"] == ["ID:A", "ID:C"] for r in union)
        assert all("path" not in r for r in union)
        assert {(r["source"], r["target"]) for r in union} == {
            ("ID:A", "ID:B"),
            ("ID:B", "ID:C"),
            ("ID:A", "ID:D"),
            ("ID:D", "ID:C"),
        }
    else:
        assert any(r.get("path") == ["ID:A", "ID:B", "ID:C"] for r in reasons)


def test_multiple_actual_path_reasons_are_retained(resources):
    extra = resources.iloc[[0]].copy()
    extra["source"], extra["target"] = "ID:B", "ID:D"
    network = Network(
        initial_nodes=["ID:A", "ID:C"], resources=pd.concat([resources, extra])
    )
    network.complete_connection(
        path_policy="all_bounded", reuse_policy="none", maxlen=3
    )
    reasons = network.list_states()[-1]["metadata"]["transition"][
        "inclusion_reasons"
    ]
    paths = [
        r["path"]
        for r in reasons
        if r["source"] == "ID:A" and r["target"] == "ID:B" and "path" in r
    ]
    assert ["ID:A", "ID:B", "ID:C"] in paths
    assert ["ID:A", "ID:B", "ID:D", "ID:C"] in paths


def test_recovery_export_and_pruned_ancestry(network, resources, tmp_path):
    network.connect_nodes()
    parent = network.current_state_id
    branch = network.edit_interactions([edit()])["state_id"]
    network.checkout(parent)
    network.edit_interactions([edit(effect="stimulation")])
    path = tmp_path / "history.json"
    network.export_history(path)
    restored = Network.from_history(path)
    assert history_payload(restored) == history_payload(network)
    assert_objects(restored)
    restored.checkout(branch)
    assert_objects(restored)
    with pytest.raises(ValueError, match="resources"):
        restored.connect_nodes()
    restored.edit_interactions([edit("set_effect", effect="bimodal")])
    online_capable = Network.from_history(path, resources=resources)
    online_capable.connect_nodes()
    network.set_max_history(2)
    assert network._pruned_state_ids
    current = network._states[network.current_state_id]
    assert current.metadata["transition"]["parent_state_id"] == parent
    assert parent in current.metadata["ancestry_gaps"]
    network.export_history(path)
    restored = Network.from_history(path)
    assert history_payload(restored) == history_payload(network)
    assert restored.history_graph().graph["pruned_state_ids"]


def test_sif_bnet_companions_and_visualization(network, resources, tmp_path):
    network.connect_nodes()
    network.edit_interactions([edit(pmids=[])])
    exporter = Exports(network)
    sif = tmp_path / "model.sif"
    exporter.export_sif(sif)
    imported = Network(sif_file=sif, resources=resources)
    manual = imported.edges.loc[
        (imported.edges.source == "ID:A") & (imported.edges.target == "ID:C")
    ].iloc[0]
    assert manual.Provenance["origin"] == "manual"
    assert manual.Provenance["edits"][0]["evidence"]["pmids"] == []
    sidecar = json.loads((tmp_path / "model.sif.neko.json").read_text())
    assert sidecar["export"]["state_id"] == network.current_state_id
    assert sidecar["history"]["database_name"] == "SIGNOR"
    assert history_payload(
        Network.from_history(tmp_path / "model.sif.neko.json")
    ) == history_payload(network)
    exporter.export_bnet(tmp_path / "model.bnet")
    bnet_sidecar = json.loads((tmp_path / "model_1.bnet.neko.json").read_text())
    assert bnet_sidecar["export"]["identifiers"][0] == {
        "identifier": "ID:A",
        "label": "A",
        "exported_identifier": "A",
    }
    assert history_payload(
        Network.from_history(tmp_path / "model_1.bnet.neko.json")
    ) == history_payload(network)
    html = network.history_html()
    assert "<svg" in html and "neko-state-" in html
    assert "Explicit literature-supported" in html
    assert "current" in network.history_digraph().source


def test_database_identity_survives_loader_and_default_changes(
    resources, monkeypatch
):
    monkeypatch.setitem(
        _universe._METHODS, "omnipath", lambda **_: resources.copy()
    )
    universe = _universe.omnipath()
    network = Network(resources=universe)
    assert network.database_name == "OmniPath"
    universe.database_name = "SIGNOR"
    assert network.database_name == "OmniPath"
    assert network.list_states()[0]["metadata"]["database_name"] == "OmniPath"
    with pytest.raises(AttributeError):
        network.database_name = "SIGNOR"
    with pytest.raises(ValueError, match="conflicts"):
        Network(resources="omnipath", database_name="SIGNOR")
    assert Network(resources=resources).database_name is None
    monkeypatch.setattr(
        _universe._signor, "signor", lambda **_: resources.copy()
    )
    assert Network(resources=_universe.signor()).database_name == "SIGNOR"
    assert (
        Network(resources=_universe.signor().interactions).database_name
        == "SIGNOR"
    )


def test_history_reload_rejects_broken_graph(network, tmp_path):
    network.edit_interactions([edit()])
    payload = history_payload(network)
    payload["states"][1]["parent_ids"] = [999]
    path = tmp_path / "broken.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="relationship"):
        Network.from_history(path)


def test_parallel_resource_evidence_and_object_effects_are_preserved(
    network, resources
):
    parallel = resources.iloc[[0]].copy()
    parallel["references"] = "999"
    parallel["is_stimulation"], parallel["is_inhibition"] = False, True
    net = Network(
        initial_nodes=["ID:A", "ID:B"],
        resources=pd.concat([resources, parallel]),
        database_name="SIGNOR",
    )
    net.connect_nodes()
    row = net.edges.iloc[0]
    assert (
        row.Effect == "stimulation"
    )  # Existing first-row policy is unchanged.
    assert {item["references"] for item in row.Provenance["originals"]} == {
        "100",
        "999",
    }
    net.edit_interactions([edit("set_effect", target="B")])
    assert {
        item["references"] for item in net.edges.iloc[0].Provenance["originals"]
    } == {"100", "999"}
    duplicate = net.edges.iloc[[0]].copy()
    duplicate["Effect"] = "form complex"
    net.edges = pd.concat([net.edges, duplicate], ignore_index=True)
    net.sync_edges_from_df()
    assert len(net.edges_as_objects()) == 2


def test_long_explanations_untracked_changes_and_working_state_roundtrip(
    network, tmp_path
):
    request = edit()
    request["rationale"] = "<script>not executable</script>" + "reason " * 100
    network.edit_interactions([request])
    html = network.history_html()
    assert "<script>not executable</script>" not in html
    assert "&lt;script&gt;not executable&lt;/script&gt;" in html
    assert (
        request["rationale"].strip()
        == network.list_states()[-1]["metadata"]["manual_changes"][0]["after"][
            "Provenance"
        ]["edits"][0]["rationale"]
    )
    with network.suspend_history():
        network.add_node("ID:D")
    path = tmp_path / "untracked.json"
    network.export_history(path)
    payload = json.loads(path.read_text())
    assert payload["working_state_matches_current"] is False
    restored = Network.from_history(path)
    assert history_payload(restored) == payload
    network.edit_interactions([edit("set_effect", effect="stimulation")])
    transition = network.list_states()[-1]["metadata"]["transition"]
    assert transition["explanation_status"] == "partial"
    assert any("already differed" in item for item in transition["limitations"])
    assert transition["nodes"]["added"][0]["Uniprot"] == "ID:D"


@pytest.mark.parametrize(
    "operation", ["radial", "cascade", "subgroup", "component"]
)
def test_other_construction_paths_are_attributed(resources, operation):
    net = Network(initial_nodes=["ID:A", "ID:C"], resources=resources)
    if operation == "radial":
        net.connect_network_radially(max_len=2, direction="OUT")
    elif operation == "cascade":
        net.connect_to_upstream_nodes(["ID:C"], depth=2, rank=2)
    elif operation == "subgroup":
        net.connect_subgroup(["ID:A", "ID:C"], maxlen=2)
    else:
        net.connect_component(["ID:A"], ["ID:C"], maxlen=2)
    transition = net.list_states()[-1]["metadata"]["transition"]
    assert transition["edges"]["added"]
    assert transition["operations"]
    for row in transition["edges"]["added"]:
        supporting = [
            reason
            for reason in transition["inclusion_reasons"]
            if (reason["source"], reason["target"])
            == (row["source"], row["target"])
        ]
        assert supporting and all(
            reason["kind"] != "unavailable" for reason in supporting
        )
        assert all(reason["operation_indices"] for reason in supporting)
    assert_objects(net)


def test_effect_change_on_legacy_row_with_empty_provenance_preserves_old_effect(
    network,
):
    network.connect_nodes()
    network.edges.at[0, "Provenance"] = None
    result = network.edit_interactions([edit("set_effect", target="B")])
    prov = result["changes"][0]["after"]["Provenance"]
    assert prov["origin"] == "unknown"
    assert prov["original"]["effect"] == "stimulation"
    assert prov["edits"][0]["before_effect"] == "stimulation"


def test_two_manual_edits_have_one_revision_and_returned_records_are_independent(
    network,
):
    parent = network.current_state_id
    result = network.edit_interactions([edit(), edit(source="C", target="A")])
    assert result["state_id"] == parent + 1
    assert len(network.list_states()) == 2
    assert len(result["changes"]) == 2
    assert all(
        change["after"]["Provenance"]["edits"][0]["state_id"]
        == result["state_id"]
        for change in result["changes"]
    )
    result["changes"][0]["after"]["Provenance"]["edits"].clear()
    assert network.edges.iloc[0].Provenance["edits"]


def test_older_serialized_network_reports_unknown_explanations(
    network, tmp_path
):
    network.connect_nodes()
    for name in (
        "_database_name",
        "_pruned_state_ids",
        "_transition_context",
        "_resources_available",
    ):
        delattr(network, name)
    for state in network._states.values():
        state.metadata.pop("transition", None)
        state.metadata.pop("database_name", None)
        state.metadata.pop("initial_nodes", None)
    restored = pickle.loads(pickle.dumps(network))
    assert restored.database_name is None
    assert all(
        state["metadata"]["transition"]["explanation_status"] == "unavailable"
        for state in restored.list_states()
    )
    restored.checkout(restored.current_state_id)
    assert_objects(restored)
    restored.export_history(tmp_path / "migrated.json")
    reloaded = Network.from_history(tmp_path / "migrated.json")
    assert history_payload(reloaded) == history_payload(restored)


def test_public_state_inspection_is_readonly_and_does_not_checkout(network):
    initial = network.current_state_id
    network.edit_interactions([edit()])
    current = network.current_state_id
    saved = network.get_state(initial)
    assert saved["id"] == initial and not saved["current"]
    assert saved["edges"]["records"] == []
    assert network.current_state_id == current
    payload = network.history_dict()
    payload["states"][-1]["metadata"]["manual_changes"].clear()
    assert network.get_state()["metadata"]["manual_changes"]
    with pytest.raises(ValueError):
        network.get_state(True)
    with pytest.raises(ValueError):
        network.get_state(999)


@pytest.mark.parametrize(
    "group",
    [
        lambda net: np.array(["ID:A", "ID:B"]),
        lambda net: net.nodes["Uniprot"].iloc[:2],
        lambda net: pd.Index(["ID:A", "ID:B"]),
        lambda net: frozenset({"ID:A", "ID:B"}),
        lambda net: (node for node in ["ID:A", "ID:B"]),
    ],
)
def test_array_like_arguments_are_recorded_not_rejected(network, group):
    network.connect_subgroup(group(network))
    assert len(network.edges) == 1
    operation = network.get_state()["metadata"]["transition"]["operations"][0]
    assert operation["method"] == "connect_subgroup"
    json.dumps(operation["parameters"], allow_nan=False)


def test_failed_history_recording_does_not_leak_transition_context(
    network, monkeypatch
):
    import neko.core.history as history

    monkeypatch.setattr(
        history,
        "bound_parameters",
        lambda *_, **__: (_ for _ in ()).throw(RuntimeError("boom")),
    )
    # add_node only records its parameters after it changed the network.
    with pytest.raises(RuntimeError):
        network.add_node("ID:D")
    assert network._transition_context is None


def test_sif_reload_does_not_turn_missing_references_into_text(
    network, resources, tmp_path
):
    network.edit_interactions([edit(pmids=[])])
    assert network.edges.iloc[0].References is None
    sif = tmp_path / "model.sif"
    Exports(network).export_sif(sif)
    assert "# Reference PMID: None" in sif.read_text()
    imported = Network(sif_file=sif, resources=resources)
    assert imported.edges.iloc[0].References == "SIF file"


def test_copy_frame_keeps_nested_values_with_duplicate_labels():
    from neko.core.history import copy_frame

    frame = pd.DataFrame(
        {"source": ["A", "B"], "Provenance": [{"x": [1]}, {"y": [2]}]},
        index=[0, 0],
    )
    copied = copy_frame(frame)
    assert copied.Provenance.tolist() == [{"x": [1]}, {"y": [2]}]
    copied.Provenance.iloc[0]["x"].append(3)
    assert frame.Provenance.iloc[0] == {"x": [1]}


def test_network_pickled_before_history_existed_can_still_grow(network):
    network.add_node("ID:D")
    for name in [
        name for name in vars(network)
        if name.startswith(("_state", "_current_state", "_root_state"))
        or name in {
            "_auto_state_depth", "_history_enabled", "_max_history",
            "_pruned_state_ids", "_transition_context",
            "_resources_available", "_database_name",
        }
    ]:
        delattr(network, name)
    restored = pickle.loads(pickle.dumps(network))
    restored.connect_nodes()
    assert len(restored.edges) == 4
    root = restored.get_state()
    assert root["id"] == restored.root_state_id
    assert root["metadata"]["transition"]["explanation_status"] == "partial"


def test_history_displays_stay_short_and_keep_user_metadata(network, capsys):
    network.connect_nodes()
    network.save_state(metadata={"note": "after manual pruning"})
    network.describe_history()
    printed = capsys.readouterr().out
    assert "State 1 - {'method': 'connect_nodes', 'args': [], 'kwargs': {}}" in printed
    assert "transition" not in printed and "Provenance" not in printed
    source = network.history_digraph().source
    assert "note: after manual pruning" in source
    assert "database_name: SIGNOR" in source
    assert "transition" not in source


def test_history_files_follow_the_umask(network, tmp_path):
    import os
    import stat

    previous = os.umask(0o022)
    try:
        network.export_history(tmp_path / "history.json")
        Exports(network).export_sif(tmp_path / "model.sif")
    finally:
        os.umask(previous)
    for name in ("history.json", "model.sif", "model.sif.neko.json"):
        assert stat.S_IMODE((tmp_path / name).stat().st_mode) == 0o644
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "history.json", "model.sif", "model.sif.neko.json",
    ]


def test_seeded_state_fingerprints_match_recomputed_ones(network):
    from neko.core.history import frame_fingerprint, state_fingerprints

    network.connect_nodes()
    network.edit_interactions([edit()])
    network.add_node("ID:D")
    for state in network._states.values():
        cached = state_fingerprints(state)
        assert cached == (
            frame_fingerprint(state.nodes), frame_fingerprint(state.edges)
        )
    assert history_payload(network)["working_state_matches_current"]


def test_single_bnet_variant_keeps_self_contained_companion(network, tmp_path):
    # MCP export reads <stem>_1.bnet.neko.json and its inline history.
    network.connect_nodes()
    Exports(network).export_bnet(tmp_path / "model.bnet")
    companion = json.loads((tmp_path / "model_1.bnet.neko.json").read_text())
    assert companion["schema_version"] == 1
    assert companion["history"] == history_payload(network)
    assert "history_file" not in companion
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "model_1.bnet", "model_1.bnet.neko.json",
    ]


def test_bnet_variants_share_one_history_file(network, tmp_path):
    network.connect_nodes()
    network.edit_interactions([edit(effect="bimodal")])
    exporter = Exports(network)
    exporter.export_bnet(tmp_path / "model.bnet")
    shared = tmp_path / "model.bnet.history.json"
    assert json.loads(shared.read_text()) == history_payload(network)
    for variant in ("model_1.bnet", "model_2.bnet"):
        companion = json.loads((tmp_path / f"{variant}.neko.json").read_text())
        assert companion["schema_version"] == 2
        assert "history" not in companion
        assert companion["history_file"]["file"] == shared.name
        assert companion["export"]["file"] == variant
        assert companion["export"]["state_id"] == network.current_state_id
        assert history_payload(
            Network.from_history(tmp_path / f"{variant}.neko.json")
        ) == history_payload(network)
    assert history_payload(Network.from_history(shared)) == history_payload(
        network
    )
    effects = {
        json.dumps(
            json.loads((tmp_path / f"model_{i}.bnet.neko.json").read_text())[
                "export"]["interactions"]
        )
        for i in (1, 2)
    }
    assert len(effects) == 2

    single = tmp_path / "single"
    exporter.export_bnet(single / "model.bnet", n=1)
    companion = json.loads((single / "model_1.bnet.neko.json").read_text())
    assert companion["schema_version"] == 1
    assert not (single / "model.bnet.history.json").exists()


def test_shared_history_reference_is_verified(network, tmp_path):
    network.edit_interactions([edit(effect="bimodal")])
    Exports(network).export_bnet(tmp_path / "model.bnet")
    companion = tmp_path / "model_1.bnet.neko.json"
    shared = tmp_path / "model.bnet.history.json"

    shared.write_text(shared.read_text().replace("SIGNOR", "SIGNOX"))
    with pytest.raises(ValueError, match="SHA-256"):
        Network.from_history(companion)
    shared.unlink()
    with pytest.raises(FileNotFoundError, match="model.bnet.history.json"):
        Network.from_history(companion)
    payload = json.loads(companion.read_text())
    for name in ("../model.bnet.history.json", "/tmp/x.json", "..", ""):
        payload["history_file"]["file"] = name
        companion.write_text(json.dumps(payload))
        with pytest.raises(ValueError, match="Malformed shared history"):
            Network.from_history(companion)
