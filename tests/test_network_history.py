import pytest
import pandas as pd
from neko.core.network import Network
import time

def test_save_and_undo():
    net = Network(initial_nodes=["TP53", "MDM2"], resources="omnipath")
    net.save_state(metadata={"description": "After init"})
    n0 = len(net.nodes)
    net.add_node("AKT1")
    net.save_state(metadata={"description": "Added AKT1"})
    n1 = len(net.nodes)
    assert n1 == n0 + 1
    net.undo()
    n2 = len(net.nodes)
    assert n2 == n0

def test_compare_states():
    net = Network(initial_nodes=["TP53", "MDM2"], resources="omnipath")
    net.save_state(metadata={"description": "After init"})
    result = net.add_node("AKT1")
    print(f"add_node('AKT1') returned: {result}")
    print("Nodes after add_node:", net.nodes)
    net.save_state(metadata={"description": "Added AKT1"})
    diff = net.compare_states(0, 1)
    print("compare_states diff:", diff)
    assert "AKT1" in diff["added_nodes"] or "P31749" in diff["added_nodes"]
    assert not diff["removed_nodes"]

def test_provenance_tracking():
    net = Network(initial_nodes=["TP53", "MDM2"], resources="omnipath")
    provenance = {"strategy": "manual_add", "params": {"source": "TP53", "target": "MDM2"}, "timestamp": time.time()}
    edge_df = pd.DataFrame({"source": ["TP53"], "target": ["MDM2"], "type": ["activation"], "references": ["PMID:12345"]})
    net.add_edge(edge_df, provenance=provenance)
    prov = net.get_edge_provenance("TP53", "MDM2")
    assert prov is not None and "manual_add" in prov
    filtered = net.filter_edges_by_provenance("manual_add")
    assert not filtered.empty
