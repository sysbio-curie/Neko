from itertools import combinations

import pandas as pd
import pytest

from neko.core.network import Network
from neko._methods.enrichment_methods import Connections


A = "P00533"
B = "P01116"
X = "P15056"


def _canonical_edges(network):
    columns = ["source", "target", "Type", "Effect", "References"]
    return (
        network.edges[columns]
        .fillna("<NA>")
        .sort_values(columns)
        .reset_index(drop=True)
    )


def _resource_rows():
    return pd.DataFrame([
        {
            "source": A,
            "target": X,
            "is_directed": True,
            "is_stimulation": True,
            "is_inhibition": False,
            "form_complex": False,
            "references": "PMID:1",
            "type": "activation",
        },
        {
            "source": A,
            "target": X,
            "is_directed": True,
            "is_stimulation": False,
            "is_inhibition": True,
            "form_complex": False,
            "references": "PMID:2",
            "type": "inhibition",
        },
        {
            "source": X,
            "target": B,
            "is_directed": True,
            "is_stimulation": False,
            "is_inhibition": False,
            "form_complex": True,
            "references": "PMID:3",
            "type": "complex",
        },
        {
            "source": B,
            "target": A,
            "is_directed": True,
            "is_stimulation": False,
            "is_inhibition": False,
            "form_complex": False,
            "references": "PMID:4",
            "type": "unknown",
        },
    ])


def test_indexed_interaction_lookup_preserves_parallel_resource_rows():
    resources = _resource_rows()
    connections = Connections(resources)

    indexed = connections.find_interactions(A, X)
    expected = resources[(resources["source"] == A) & (resources["target"] == X)]

    pd.testing.assert_frame_equal(indexed, expected)
    assert connections.find_interactions(A, B).empty
    assert connections.find_target_neighbours(A) == [X]


@pytest.mark.parametrize("only_signed", [False, True])
def test_batched_connect_nodes_matches_legacy_edge_mutation(only_signed):
    resources = _resource_rows()
    base = Network(initial_nodes=[A, B, X], resources=resources)
    base._history_enabled = False
    legacy = base.copy()
    batched = base.copy()

    for node1, node2 in combinations(legacy.nodes["Uniprot"], 2):
        for source, target in ((node1, node2), (node2, node1)):
            interaction = legacy.resources.loc[
                (legacy.resources["source"] == source)
                & (legacy.resources["target"] == target)
            ]
            if not interaction.empty and (
                    not only_signed
                    or legacy.check_sign(interaction) != "undefined"
                ):
                legacy.add_edge(interaction)

    batched.connect_nodes(only_signed=only_signed)

    pd.testing.assert_frame_equal(
        _canonical_edges(batched),
        _canonical_edges(legacy),
    )
    assert {
        (edge.source, edge.target, edge.interaction_type, edge.metadata["Effect"])
        for edge in batched.edges_as_objects()
    } == {
        (edge.source, edge.target, edge.interaction_type, edge.metadata["Effect"])
        for edge in legacy.edges_as_objects()
    }


def test_unsigned_working_edge_is_not_a_signed_path():
    working_edges = pd.DataFrame([{
        "source": A,
        "target": B,
        "Type": "interaction",
        "Effect": "undefined",
        "References": "PMID:1",
    }])

    connections = Connections(working_edges)

    assert connections.bfs(
        A,
        B,
        maxlen=1,
        only_signed=True,
    ) == []
