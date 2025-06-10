import pytest
import pandas as pd
from neko.core.network import Network
import os
import difflib
from neko._outputs.exports import Exports

@pytest.fixture
def sample_genes():
    # Use Uniprot IDs for test genes
    return ["P12931", "P46531", "P12830", "P19022"]  # SRC, NOTCH1, CDH1, CDH2

@pytest.fixture
def sample_resources():
    # Minimal mock resource DataFrame with Uniprot IDs
    return pd.DataFrame({
        "source": ["P12931", "P46531", "P12830"],
        "target": ["P12830", "P19022", "P19022"],
        "Type": ["activation", "inhibition", "activation"],
        "Effect": ["stimulation", "inhibition", "stimulation"],
        "References": ["PMID:1", "PMID:2", "PMID:3"]
    })

def test_network_creation_from_genes(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    assert set(net.nodes["Uniprot"]).intersection(sample_genes)
    assert net.resources is not None

def test_add_and_remove_node(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    # Add a node that is present in resources
    present_uniprot = "P12931"  # already in resources
    net.add_node(present_uniprot)
    assert present_uniprot in net.nodes["Uniprot"].values
    net.remove_node(present_uniprot)
    assert present_uniprot not in net.nodes["Uniprot"].values

    # Try to add a node not present in resources
    absent_uniprot = "Q9Y2X3"
    net.add_node(absent_uniprot)
    assert absent_uniprot not in net.nodes["Uniprot"].values

def test_add_and_remove_edge(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    edge_df = pd.DataFrame({
        "source": ["P12931"],
        "target": ["P19022"],
        "type": ["activation"],
        "references": ["PMID:123"]
    })
    net.add_edge(edge_df)
    assert ((net.edges["source"] == "P12931") & (net.edges["target"] == "P19022")).any()
    net.remove_edge("P12931", "P19022")
    assert not ((net.edges["source"] == "P12931") & (net.edges["target"] == "P19022")).any()

def test_connect_nodes_and_complete_connection(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    net.connect_nodes(only_signed=True, consensus_only=False)
    net.complete_connection(maxlen=2, algorithm="dfs", only_signed=True)
    assert isinstance(net.edges, pd.DataFrame)

def test_remove_bimodal_and_undefined(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    # Add bimodal and undefined edges
    net.edges.loc[len(net.edges)] = ["P12931", "P19022", "type", "bimodal", "PMID:4"]
    net.edges.loc[len(net.edges)] = ["P12931", "P12830", "type", "undefined", "PMID:5"]
    net.remove_bimodal_interactions()
    assert not (net.edges["Effect"] == "bimodal").any()
    net.remove_undefined_interactions()
    assert not (net.edges["Effect"] == "undefined").any()

def test_check_nodes_and_connectivity(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    present = net.check_nodes(["P12931", "P12830", "Q99999"])
    assert "P12931" in present and "Q99999" not in present
    assert isinstance(net.is_connected(), bool)

def test_network_omnipath_benchmark(tmp_path):
    # 1. Define the gene list
    genes = ["SRC", "NOTCH1", "FAK", "CDH1", "CDH2", "VIM"]

    # 2. Create the network with Omnipath as resource
    net = Network(initial_nodes=genes, resources="omnipath")

    # 3. Complete connection (fully connect the network) with explicit parameters
    net.complete_connection(maxlen=2, algorithm="bfs", minimal=True, only_signed=True, consensus=False, connect_with_bias=False)

    # 4. Export to SIF (temporary output)
    out_sif = tmp_path / "benchmark_network.sif"
    Exports(net).export_sif(str(out_sif))

    # 5. Compare with reference SIF if it exists, otherwise create it
    ref_sif = os.path.join(os.path.dirname(__file__), "benchmark_network.sif")
    if os.path.exists(ref_sif):
        with open(ref_sif) as f1, open(out_sif) as f2:
            ref_lines = f1.readlines()
            out_lines = f2.readlines()
        if ref_lines != out_lines:
            diff = list(difflib.unified_diff(ref_lines, out_lines, fromfile="reference", tofile="current"))
            print("\n".join(diff))
            assert False, "benchmark_network.sif differs from reference. See diff above."
    else:
        # Save the new SIF as reference for future runs
        import shutil
        shutil.copy(str(out_sif), ref_sif)
        print(f"Reference SIF file created at {ref_sif}. Please review and commit if correct.")
