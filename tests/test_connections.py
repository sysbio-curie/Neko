import pandas as pd
import pytest

from neko.core.network import Network
from neko.core.strategies import complete_connection
from neko.core.strategy_options import (
    ConnectionStrategyMigrationWarning,
    resolve_connection_policies,
)
from neko._methods.enrichment_methods import Connections


A = "P00533"  # EGFR
B = "P01116"  # KRAS
X = "P15056"  # BRAF
Y = "Q02750"  # MAP2K1
Z = "P28482"  # MAPK1
W = "P01100"  # FOS


def _resource(edges):
    """Build a small signed resource with stable, distinguishable evidence."""

    rows = []
    for index, (source, target) in enumerate(edges, start=1):
        rows.append({
            "source": source,
            "target": target,
            "is_directed": True,
            "is_stimulation": True,
            "is_inhibition": False,
            "form_complex": False,
            "references": f"PMID:{index}",
            "type": "activation",
        })
    return pd.DataFrame(rows)


@pytest.fixture
def alternative_path_resource():
    return _resource([
        (A, X),
        (X, B),
        (A, Y),
        (Y, B),
        (A, Z),
        (Z, W),
        (W, B),
    ])


def _edge_keys(network):
    return set(network.edges[["source", "target"]].itertuples(
        index=False,
        name=None,
    ))


def _canonical_edges(network):
    columns = ["source", "target", "Type", "Effect", "References"]
    return (
        network.edges[columns]
        .fillna("<NA>")
        .sort_values(columns)
        .reset_index(drop=True)
    )


def test_bfs_returns_one_shortest_path_with_an_inclusive_cutoff(
        alternative_path_resource):
    connections = Connections(alternative_path_resource)

    paths = connections.bfs(A, B, maxlen=2, only_signed=True)

    assert len(paths) == 1
    assert tuple(paths[0]) in {
        (A, X, B),
        (A, Y, B),
    }
    assert connections.bfs(A, B, maxlen=1, only_signed=True) == []


def test_dfs_returns_all_simple_paths_through_the_inclusive_cutoff(
        alternative_path_resource):
    connections = Connections(alternative_path_resource)

    paths = connections.find_paths(
        A,
        B,
        maxlen=3,
        only_signed=True,
    )

    assert {tuple(path) for path in paths} == {
        (A, X, B),
        (A, Y, B),
        (A, Z, W, B),
    }


def test_all_shortest_returns_the_shortest_path_edge_union(
        alternative_path_resource):
    connections = Connections(alternative_path_resource)

    edges = connections.bfs_all_shortest_edges(
        A,
        B,
        maxlen=3,
        only_signed=True,
    )

    assert set(edges) == {
        (A, X),
        (X, B),
        (A, Y),
        (Y, B),
    }
    assert (A, Z) not in edges


def test_connect_nodes_builds_the_induced_subgraph_without_bridge_nodes(
        alternative_path_resource):
    network = Network(
        initial_nodes=[A, B, X, Y],
        resources=alternative_path_resource,
    )
    network._history_enabled = False

    network.connect_nodes(only_signed=True)

    assert set(network.nodes["Uniprot"]) == {A, B, X, Y}
    assert _edge_keys(network) == {
        (A, X),
        (X, B),
        (A, Y),
        (Y, B),
    }


def test_legacy_bfs_completion_adds_one_shortest_route(
        alternative_path_resource):
    network = Network(
        initial_nodes=[A, B],
        resources=alternative_path_resource,
    )
    network._history_enabled = False

    with pytest.warns(ConnectionStrategyMigrationWarning):
        complete_connection(
            network,
            maxlen=3,
            algorithm="bfs",
            minimal=True,
            only_signed=True,
            connect_with_bias=False,
        )

    edges = _edge_keys(network)
    assert edges in (
        {(A, X), (X, B)},
        {(A, Y), (Y, B)},
    )


def test_legacy_dfs_completion_adds_every_bounded_route(
        alternative_path_resource):
    network = Network(
        initial_nodes=[A, B],
        resources=alternative_path_resource,
    )
    network._history_enabled = False

    with pytest.warns(ConnectionStrategyMigrationWarning):
        complete_connection(
            network,
            maxlen=3,
            algorithm="dfs",
            minimal=True,
            only_signed=True,
            connect_with_bias=False,
        )

    assert _edge_keys(network) == {
        (A, X),
        (X, B),
        (A, Y),
        (Y, B),
        (A, Z),
        (Z, W),
        (W, B),
    }


def test_all_shortest_completion_adds_both_shortest_routes_only(
        alternative_path_resource):
    network = Network(
        initial_nodes=[A, B],
        resources=alternative_path_resource,
    )
    network._history_enabled = False

    complete_connection(
        network,
        maxlen=3,
        path_policy="all_shortest",
        reuse_policy="discovered_paths",
        only_signed=True,
    )

    assert _edge_keys(network) == {
        (A, X),
        (X, B),
        (A, Y),
        (Y, B),
    }


def test_induced_reuse_can_avoid_an_alternative_bridge():
    resources = _resource([
        (B, Y),
        (Y, A),
        (A, Y),
        (Y, B),
        (A, X),
        (X, B),
    ])
    discovered = Network(initial_nodes=[A, B], resources=resources)
    discovered._history_enabled = False
    induced = discovered.copy()

    complete_connection(
        discovered,
        maxlen=2,
        path_policy="one_shortest",
        reuse_policy="discovered_paths",
        only_signed=True,
    )
    complete_connection(
        induced,
        maxlen=2,
        path_policy="one_shortest",
        reuse_policy="induced_subgraph",
        only_signed=True,
    )

    assert X in set(discovered.nodes["Uniprot"])
    assert X not in set(induced.nodes["Uniprot"])
    assert set(induced.nodes["Uniprot"]) == {A, B, Y}
    assert _edge_keys(induced) == {
        (B, Y),
        (Y, A),
        (A, Y),
        (Y, B),
    }


def test_discovered_path_reuse_reduces_resource_searches():
    resources = _resource([
        (A, X),
        (X, B),
    ])

    def run(reuse_policy):
        network = Network(initial_nodes=[A, B, X], resources=resources)
        network._history_enabled = False
        resource_bfs = network._connect.bfs
        calls = []

        def counted_bfs(*args, **kwargs):
            calls.append((kwargs.get("start"), kwargs.get("end")))
            return resource_bfs(*args, **kwargs)

        network._connect.bfs = counted_bfs
        complete_connection(
            network,
            maxlen=2,
            path_policy="one_shortest",
            reuse_policy=reuse_policy,
            only_signed=True,
        )
        return network, calls

    independent, independent_calls = run("none")
    reused, reused_calls = run("discovered_paths")

    assert _edge_keys(independent) == _edge_keys(reused)
    assert len(reused_calls) < len(independent_calls)


@pytest.mark.parametrize(
    ("minimal", "bias", "expected_reuse"),
    [
        (False, False, "none"),
        (True, False, "discovered_paths"),
        (False, True, "induced_subgraph"),
        (True, True, "induced_subgraph"),
    ],
)
@pytest.mark.parametrize(
    ("algorithm", "expected_path"),
    [
        ("bfs", "one_shortest"),
        ("dfs", "all_bounded"),
    ],
)
def test_legacy_reuse_arguments_map_to_explicit_policies(
        minimal, bias, expected_reuse, algorithm, expected_path):
    with pytest.warns(ConnectionStrategyMigrationWarning):
        resolved = resolve_connection_policies(
            maxlen=2,
            path_policy=None,
            reuse_policy=None,
            algorithm=algorithm,
            minimal=minimal,
            connect_with_bias=bias,
        )

    assert resolved.path_policy == expected_path
    assert resolved.reuse_policy == expected_reuse
    assert resolved.maxlen == 2


@pytest.mark.parametrize(
    ("minimal", "bias", "reuse_policy"),
    [
        (False, False, "none"),
        (True, False, "discovered_paths"),
        (False, True, "induced_subgraph"),
        (True, True, "induced_subgraph"),
    ],
)
@pytest.mark.parametrize(
    ("algorithm", "path_policy"),
    [
        ("bfs", "one_shortest"),
        ("dfs", "all_bounded"),
    ],
)
def test_legacy_and_explicit_policy_calls_build_the_same_topology(
        minimal, bias, reuse_policy, algorithm, path_policy):
    resources = _resource([
        (B, Y),
        (Y, A),
        (A, Y),
        (Y, B),
        (A, X),
        (X, B),
    ])
    legacy = Network(initial_nodes=[A, B], resources=resources)
    legacy._history_enabled = False
    explicit = legacy.copy()

    with pytest.warns(ConnectionStrategyMigrationWarning):
        complete_connection(
            legacy,
            maxlen=2,
            algorithm=algorithm,
            minimal=minimal,
            connect_with_bias=bias,
            only_signed=True,
        )
    complete_connection(
        explicit,
        maxlen=2,
        path_policy=path_policy,
        reuse_policy=reuse_policy,
        only_signed=True,
    )

    assert set(legacy.nodes["Uniprot"]) == set(explicit.nodes["Uniprot"])
    pd.testing.assert_frame_equal(
        _canonical_edges(legacy),
        _canonical_edges(explicit),
    )


def test_new_policies_require_a_finite_positive_cutoff():
    with pytest.raises(ValueError, match="positive integer"):
        resolve_connection_policies(
            maxlen=None,
            path_policy="one_shortest",
            reuse_policy="discovered_paths",
        )


def test_legacy_unbounded_bfs_maps_to_the_historical_cutoff():
    with pytest.warns(
            ConnectionStrategyMigrationWarning,
            match="maxlen=10",
        ):
        resolved = resolve_connection_policies(
            maxlen=None,
            path_policy=None,
            reuse_policy=None,
            algorithm="bfs",
        )

    assert resolved.maxlen == 10


def test_mixed_legacy_and_new_policies_are_rejected():
    with pytest.raises(ValueError, match="Do not mix"):
        resolve_connection_policies(
            maxlen=2,
            path_policy="all_shortest",
            reuse_policy="discovered_paths",
            algorithm="bfs",
        )


def test_history_does_not_translate_policy_strings():
    network = Network(
        initial_nodes=[A, B],
        resources=_resource([(A, B)]),
    )

    def unexpected_lookup(identifier):
        raise AssertionError(f"unexpected identifier lookup: {identifier}")

    network.mapping_node_identifier = unexpected_lookup
    network.complete_connection(
        maxlen=1,
        path_policy="one_shortest",
        reuse_policy="discovered_paths",
        only_signed=True,
    )

    metadata = network.list_states()[-1]["metadata"]
    assert metadata["kwargs"]["path_policy"] == "one_shortest"
    assert metadata["kwargs"]["reuse_policy"] == "discovered_paths"
