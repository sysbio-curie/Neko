import pytest
import pandas as pd
from neko.core.network import Network
from neko.core import strategies

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

def test_connect_genes_to_phenotype_omnipath():
    genes = ["SRC", "NOTCH1", "FAK"]
    net = Network(genes, resources="omnipath")
    strategies.connect_genes_to_phenotype(net, id_accession="GO:0001837", phenotype="epithelial to mesenchymal transition", only_signed=True, compress=True, maxlen=1)
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
