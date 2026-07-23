import os

import pandas as pd
import pytest
from neko.core.network import Network
from neko.core import strategies
from neko._annotations.gene_ontology import GOGene, GOTerm


class _Ontology:
    accession_to_phenotype_dict = {}

    def resolve_accession(self, phenotype=None, id_accession=None):
        return id_accession or 'GO:0062043'

    def get_term(self, go_id):
        return GOTerm(
            go_id='GO:0062043',
            label=(
                'positive regulation of cardiac epithelial to mesenchymal '
                'transition'
            ),
        )

    def fetch_go_genes(self, go_id, **kwargs):
        self.fetch_kwargs = kwargs
        return [
            GOGene(
                gene_id='UniProtKB:Q04771',
                symbol='ACVR1',
                taxon_id='NCBITaxon:9606',
            ),
            GOGene(
                gene_id='Ensembl:ENSG000001234',
                symbol='FALLBACK',
                taxon_id='NCBITaxon:9606',
            ),
        ]


class _PhenotypeNetwork:
    def __init__(self):
        self._ontology = _Ontology()
        self.nodes = pd.DataFrame([{
            'Genesymbol': 'SRC',
            'Uniprot': 'P12931',
            'Type': 'NaN',
        }])
        self.edges = pd.DataFrame(
            columns=['source', 'target', 'Type', 'Effect', 'References'],
        )
        self.mapping_calls = []

    def mapping_node_identifier(self, node):
        self.mapping_calls.append(node)
        return [None, node, 'P99999']


def test_connect_genes_to_phenotype_prefers_go_uniprot_ids(monkeypatch):
    network = _PhenotypeNetwork()
    calls = []

    def connect_component(network, comp_a, comp_b, **kwargs):
        calls.append((comp_a, comp_b, kwargs))

    monkeypatch.setattr(strategies, 'connect_component', connect_component)

    strategies.connect_genes_to_phenotype(
        network,
        id_accession='GO:0062043',
        taxon_id='NCBITaxon:9606',
        include_descendants=True,
        exclude_automatic_assertions=True,
    )

    assert network.mapping_calls == ['FALLBACK']
    assert calls[0][0] == ['P12931']
    assert set(calls[0][1]) == {'Q04771', 'P99999'}
    assert network._ontology.fetch_kwargs == {
        'taxon_id': 'NCBITaxon:9606',
        'include_descendants': True,
        'exclude_automatic_assertions': True,
    }


def test_connect_genes_to_phenotype_uses_canonical_label(monkeypatch):
    network = _PhenotypeNetwork()

    def connect_component(network, comp_a, comp_b, **kwargs):
        network.nodes.loc[len(network.nodes)] = {
            'Genesymbol': 'ACVR1',
            'Uniprot': 'Q04771',
            'Type': 'NaN',
        }
        network.edges.loc[len(network.edges)] = {
            'source': 'P12931',
            'target': 'Q04771',
            'Type': 'interaction',
            'Effect': 'stimulation',
            'References': 'PMID:1',
        }

    monkeypatch.setattr(strategies, 'connect_component', connect_component)

    strategies.connect_genes_to_phenotype(
        network,
        id_accession='GO:0062043',
        phenotype='caller supplied label',
        compress=True,
    )

    canonical = (
        'positive_regulation_of_cardiac_epithelial_to_mesenchymal_transition'
    )
    assert canonical in network.nodes['Uniprot'].values
    assert canonical in network.edges['target'].values


def test_phenotype_compression_merges_opposite_effects(monkeypatch):
    network = _PhenotypeNetwork()

    def fetch_go_genes(go_id, **kwargs):
        return [
            GOGene(
                gene_id='UniProtKB:P12931',
                symbol='SRC',
                taxon_id='NCBITaxon:9606',
            ),
            GOGene(
                gene_id='UniProtKB:Q04771',
                symbol='ACVR1',
                taxon_id='NCBITaxon:9606',
            ),
        ]

    network._ontology.fetch_go_genes = fetch_go_genes

    def connect_component(network, comp_a, comp_b, **kwargs):
        network.nodes.loc[len(network.nodes)] = {
            'Genesymbol': 'ACVR1',
            'Uniprot': 'Q04771',
            'Type': 'NaN',
        }
        network.edges.loc[len(network.edges)] = {
            'source': 'P12931',
            'target': 'Q04771',
            'Type': 'interaction',
            'Effect': 'inhibition',
            'References': 'PMID:1',
        }

    monkeypatch.setattr(strategies, 'connect_component', connect_component)

    strategies.connect_genes_to_phenotype(
        network,
        id_accession='GO:0062043',
        compress=True,
    )

    phenotype = (
        'positive_regulation_of_cardiac_epithelial_to_mesenchymal_transition'
    )
    phenotype_nodes = network.nodes[
        network.nodes['Genesymbol'] == phenotype
    ]
    phenotype_edges = network.edges[
        (network.edges['source'] == 'P12931')
        & (network.edges['target'] == phenotype)
    ]

    assert len(phenotype_nodes) == 1
    assert len(phenotype_edges) == 1
    assert phenotype_edges.iloc[0]['Effect'] == 'bimodal'
    assert phenotype_edges.iloc[0]['References'] == (
        'PMID:1; Gene Ontology: GO:0062043'
    )
    assert 'Q04771' not in network.nodes['Uniprot'].values

def test_connect_nodes_omnipath():
    # Use a small set of genes known to be in Omnipath
    genes = ["TP53", "MDM2", "CDKN1A", "RB1"]
    net = Network(genes, resources="omnipath")
    strategies.connect_nodes(net)
    assert not net.edges.empty

def test_connect_subgroup_omnipath():
    genes = ["TP53", "MDM2", "CDKN1A"]
    net = Network(genes, resources="omnipath")
    strategies.connect_subgroup(net, genes, maxlen=3, only_signed=True, consensus= True)
    assert not net.edges.empty

def test_connect_component_omnipath():
    genes = ["TP53", "MDM2", "CDKN1A", "RB1"]
    net = Network(genes, resources="omnipath")
    # Split genes into two components
    comp_A = ["TP53", "MDM2"]
    comp_B = ["CDKN1A", "RB1"]
    strategies.connect_component(net, comp_A, comp_B, maxlen=3, mode="OUT", only_signed=True, consensus=True)
    assert not net.edges.empty

def test_connect_to_upstream_nodes_omnipath():
    genes = ["TP53", "MDM2", "CDKN1A", "RB1"]
    net = Network(genes, resources="omnipath")
    strategies.connect_to_upstream_nodes(net, depth=2, rank=2, only_signed=True, consensus=True)
    assert not net.edges.empty

@pytest.mark.live
@pytest.mark.skipif(
    os.environ.get('NEKO_RUN_LIVE_TESTS') != '1',
    reason='set NEKO_RUN_LIVE_TESTS=1 to call GO and OmniPath services',
)
def test_connect_genes_to_phenotype_omnipath():
    genes = ["SRC", "NOTCH1", "FAK"]
    net = Network(genes, resources="omnipath")
    strategies.connect_genes_to_phenotype(net, id_accession="GO:0062043", only_signed=True, compress=True, maxlen=1)
    assert not net.edges.empty

def test_connect_network_radially_omnipath():
    genes = ["TP53", "MDM2", "CDKN1A", "RB1"]
    net = Network(genes, resources="omnipath")
    strategies.connect_network_radially(net, max_len=1, direction=None, loops=False, consensus=True, only_signed=True)
    assert not net.edges.empty

def test_connect_as_atopo_omnipath():
    genes = ["SRC", "NOTCH1", "FAK"]
    net = Network(genes, resources="omnipath")
    strategies.connect_as_atopo(net, strategy="radial", max_len=1, loops=False, outputs=["AKT1"], only_signed=True, consensus=True)
    assert not net.edges.empty

# Add more tests for each strategy as you modularize them
