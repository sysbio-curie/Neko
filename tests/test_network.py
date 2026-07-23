import pytest
import pandas as pd

from neko.core.network import Network
import os
import difflib
from neko._outputs.exports import Exports
from neko.core.tools import is_connected

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
        "References": ["PMID:1", "PMID:2", "PMID:3"],
    })

def test_network_creation_from_genes(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    assert set(net.nodes["Uniprot"]).intersection(sample_genes)
    assert net.resources is not None


def test_aliases_share_one_canonical_visualizer_node(monkeypatch):
    import neko.core.network as network_module
    from neko._visual.visualize_network import NetworkVisualizer

    identifiers = {
        'FAK': [None, 'PTK2', 'Q05397'],
        'PTK2': [None, 'PTK2', 'Q05397'],
        'Q05397': [None, 'PTK2', 'Q05397'],
        'SRC': [None, 'SRC', 'P12931'],
        'P12931': [None, 'SRC', 'P12931'],
    }
    monkeypatch.setattr(
        network_module,
        'mapping_node_identifier',
        lambda node: identifiers[node],
    )
    resources = pd.DataFrame({
        'source': ['Q05397'],
        'target': ['P12931'],
        'Type': ['activation'],
        'Effect': ['stimulation'],
        'References': ['PMID:1'],
    })

    net = Network(
        initial_nodes=['FAK', 'PTK2', 'SRC'],
        resources=resources,
    )

    assert net.initial_nodes == ['PTK2', 'SRC']
    assert net.nodes['Genesymbol'].tolist() == ['PTK2', 'SRC']

    visualizer = NetworkVisualizer(net, noi=True)
    visualizer.tissue_mapping(pd.DataFrame({
        'Genesymbol': ['PTK2', 'SRC'],
        'in_tissue': [True, True],
    }))
    visualizer._NetworkVisualizer__build_graph()

    assert visualizer.graph.source.count('\n\tPTK2 [') == 1
    assert visualizer.graph.source.count('\n\tSRC [') == 1

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
        "references": ["PMID:123"],
    })
    net.add_edge(edge_df)
    assert ((net.edges["source"] == "P12931") & (net.edges["target"] == "P19022")).any()
    net.remove_edge("P12931", "P19022")
    assert not ((net.edges["source"] == "P12931") & (net.edges["target"] == "P19022")).any()


def test_add_edge_merges_opposite_signs_and_evidence(
        sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    stimulation = pd.DataFrame({
        "source": ["P12931"],
        "target": ["P19022"],
        "type": ["activation"],
        "references": ["PMID:stim"],
        "is_stimulation": [True],
        "is_inhibition": [False],
    })
    inhibition = pd.DataFrame({
        "source": ["P12931"],
        "target": ["P19022"],
        "type": ["inhibition"],
        "references": ["PMID:inhib"],
        "is_stimulation": [False],
        "is_inhibition": [True],
    })

    net.add_edge(stimulation)
    net.add_edge(inhibition)

    matching = net.edges[
        (net.edges["source"] == "P12931")
        & (net.edges["target"] == "P19022")
    ]
    assert len(matching) == 1
    assert matching.iloc[0]["Effect"] == "bimodal"
    assert matching.iloc[0]["Type"] == "activation; inhibition"
    assert matching.iloc[0]["References"] == "PMID:stim; PMID:inhib"
    assert len(net.edges_as_objects()) == len(net.edges)


def test_gene_symbol_conversion_prefers_custom_network_nodes(monkeypatch):
    import neko.core.network as network_module

    phenotype = "cell_cycle_arrest"
    net = Network.__new__(Network)
    net.nodes = pd.DataFrame([
        {"Genesymbol": "A", "Uniprot": "UP_A", "Type": "NaN"},
        {
            "Genesymbol": phenotype,
            "Uniprot": phenotype,
            "Type": "phenotype",
        },
    ])
    net.edges = pd.DataFrame([{
        "source": "UP_A",
        "target": phenotype,
        "Type": "interaction",
        "Effect": "stimulation",
        "References": "PMID:1",
    }])

    monkeypatch.setattr(
        network_module,
        "mapping_node_identifier",
        lambda identifier: pytest.fail(
            f"unexpected external translation for {identifier}",
        ),
    )

    converted = net.convert_edgelist_into_genesymbol()

    assert converted.loc[0, "source"] == "A"
    assert converted.loc[0, "target"] == phenotype
    assert converted[["source", "target"]].notna().all().all()


def test_gene_symbol_conversion_preserves_unknown_identifier(monkeypatch):
    import neko.core.network as network_module

    net = Network.__new__(Network)
    net.nodes = pd.DataFrame([{
        "Genesymbol": "A",
        "Uniprot": "UP_A",
        "Type": "NaN",
    }])
    net.edges = pd.DataFrame([{
        "source": "UP_A",
        "target": "custom_target",
        "Type": "interaction",
        "Effect": "stimulation",
        "References": "PMID:1",
    }])
    monkeypatch.setattr(
        network_module,
        "mapping_node_identifier",
        lambda identifier: [None, None, None],
    )

    converted = net.convert_edgelist_into_genesymbol()

    assert converted.loc[0, "source"] == "A"
    assert converted.loc[0, "target"] == "custom_target"


def test_sif_import_merges_opposite_signs(tmp_path, monkeypatch):
    import neko.core.network as network_module

    monkeypatch.setattr(
        network_module,
        "mapping_node_identifier",
        lambda identifier: [None, identifier, identifier],
    )
    monkeypatch.setattr(
        network_module,
        "check_gene_list_format",
        lambda identifiers: True,
    )
    sif_file = tmp_path / "opposite.sif"
    sif_file.write_text("A stimulation B\nA inhibition B\n")
    resources = pd.DataFrame(columns=[
        "source",
        "target",
        "is_directed",
        "is_stimulation",
        "is_inhibition",
        "form_complex",
    ])

    net = Network(sif_file=str(sif_file), resources=resources)

    assert len(net.edges) == 1
    assert net.edges.loc[0, "Effect"] == "bimodal"
    assert net.edges.loc[0, "References"] == "SIF file"

def test_connect_nodes_and_complete_connection(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    net.connect_nodes(only_signed=True, consensus_only=False)
    net.complete_connection(maxlen=2, algorithm="dfs", only_signed=True)
    assert isinstance(net.edges, pd.DataFrame)

def test_remove_bimodal_and_undefined(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    # Add bimodal and undefined edges
    net.edges.loc[len(net.edges)] = {"source": "P12931", "target": "P19022", "Type": "type", "Effect": "bimodal", "References": "PMID:4"}
    net.edges.loc[len(net.edges)] = {"source": "P12931", "target": "P12830", "Type": "type", "Effect": "undefined", "References": "PMID:5"}
    net.remove_bimodal_interactions()
    assert not (net.edges["Effect"] == "bimodal").any()
    net.remove_undefined_interactions()
    assert not (net.edges["Effect"] == "undefined").any()

def test_check_nodes_and_connectivity(sample_genes, sample_resources):
    net = Network(initial_nodes=sample_genes, resources=sample_resources)
    present = net.check_nodes(["P12931", "P12830", "Q99999"])
    assert "P12931" in present and "Q99999" not in present
    assert isinstance(is_connected(net), bool)
