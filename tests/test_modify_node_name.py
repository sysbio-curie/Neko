import pandas as pd
import pytest

import neko.core.network as network_module
from neko.core.network import Network


@pytest.fixture
def network(monkeypatch):
    identifiers = {
        "OLD": [None, "OLD", "P00001"],
        "P00001": [None, "OLD", "P00001"],
        "Q00002": [None, "OLD", "Q00002"],
        "SRC": [None, "SRC", "P12931"],
        "P12931": [None, "SRC", "P12931"],
    }
    monkeypatch.setattr(
        network_module,
        "mapping_node_identifier",
        lambda node: identifiers.get(node, [None, node, node]),
    )
    resources = pd.DataFrame({
        "source": ["P00001", "Q00002"],
        "target": ["P12931", "P12931"],
        "Type": ["activation", "activation"],
        "Effect": ["stimulation", "stimulation"],
        "References": ["PMID:1", "PMID:2"],
    })
    network = Network(
        initial_nodes=["P00001", "Q00002", "P12931"],
        resources=resources,
    )
    network._add_resource_interactions(resources)
    return network


def test_modify_node_name_changes_labels_only(network):
    identifiers_before = network.nodes["Uniprot"].copy()
    edges_before = network.edges.copy(deep=True)
    initial_nodes_before = network.initial_nodes.copy()

    network.modify_node_name("OLD", "CUSTOM")

    renamed = network.nodes[network.nodes["Genesymbol"] == "CUSTOM"]
    assert set(renamed["Uniprot"]) == {"P00001", "Q00002"}
    pd.testing.assert_series_equal(
        network.nodes["Uniprot"],
        identifiers_before,
    )
    pd.testing.assert_frame_equal(network.edges, edges_before)
    assert network.initial_nodes == initial_nodes_before
    assert {
        node.metadata["Genesymbol"]
        for node in network.nodes_as_objects()
        if node.id in {"P00001", "Q00002"}
    } == {"CUSTOM"}


@pytest.mark.parametrize("mode", ["Uniprot", "both"])
def test_modify_node_name_rejects_identifier_changes(network, mode):
    nodes_before = network.nodes.copy(deep=True)
    edges_before = network.edges.copy(deep=True)

    with pytest.raises(ValueError, match="identifiers cannot be renamed"):
        network.modify_node_name("OLD", "CUSTOM", type=mode)

    pd.testing.assert_frame_equal(network.nodes, nodes_before)
    pd.testing.assert_frame_equal(network.edges, edges_before)


def test_modify_node_name_rejects_label_collision(network):
    with pytest.raises(ValueError, match="already used"):
        network.modify_node_name("OLD", "SRC")


@pytest.mark.parametrize(
    ("old_name", "new_name", "message"),
    [
        ("MISSING", "CUSTOM", "is not in the network"),
        ("OLD", "", "new_name must be a non-empty string"),
        ("OLD", " CUSTOM", "leading or trailing whitespace"),
    ],
)
def test_modify_node_name_rejects_invalid_labels(
        network, old_name, new_name, message):
    with pytest.raises(ValueError, match=message):
        network.modify_node_name(old_name, new_name)
