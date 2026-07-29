import pandas as pd

from neko.core import tools
from neko.core.network import Network
from neko.inputs._db.psp import psp


def phosphosite_frames():
    kinase_substrate = pd.DataFrame({
        'KIN_ORGANISM': ['human', 'human'],
        'SUB_ORGANISM': ['human', 'human'],
        'GENE': ['K1', 'S1'],
        'SUB_GENE': ['S1', 'S2'],
        'SUB_MOD_RSD': ['S9', 'T10'],
    })
    regulatory_sites = pd.DataFrame({
        'ORGANISM': ['human', 'human'],
        'GENE': ['S1', 'S2'],
        'MOD_RSD': ['S9-p', 'T10-p'],
        'ON_FUNCTION': ['activity induced', 'activity induced'],
    })

    return kinase_substrate, regulatory_sites


def test_phosphosite_identifier_bypasses_uniprot(monkeypatch):
    monkeypatch.setattr(
        tools.mapping,
        'to_uniprot',
        lambda value: (_ for _ in ()).throw(
            AssertionError('phosphosites must not reach UniProt'),
        ),
    )

    assert tools.mapping_node_identifier('MAP3K4_T1494') == [
        None,
        'MAP3K4_T1494',
        'MAP3K4_T1494',
    ]
    assert tools.mapping_node_identifier('MAP3K4_t1494') == [
        None,
        'MAP3K4_T1494',
        'MAP3K4_T1494',
    ]
    assert tools.mapping_node_identifier('GENE_PART_s42') == [
        None,
        'GENE_PART_S42',
        'GENE_PART_S42',
    ]
    assert tools.check_gene_list_format(['MAP3K4_T1494']) is True


def test_phosphosite_network_uses_resource_identifiers(monkeypatch):
    kinase_substrate, regulatory_sites = phosphosite_frames()
    resources = psp(
        'human',
        kinase_substrate=kinase_substrate,
        regulatory_sites=regulatory_sites,
    )
    calls = []

    def to_uniprot(value):
        calls.append(value)
        return {'K1': 'UP_K1', 'S1': 'UP_S1', 'S2': 'UP_S2'}.get(value)

    monkeypatch.setattr(tools.mapping, 'to_uniprot', to_uniprot)
    monkeypatch.setattr(
        tools.mapping,
        'to_genesymbol',
        lambda value: {
            'UP_K1': 'K1',
            'UP_S1': 'S1',
            'UP_S2': 'S2',
        }.get(value, value),
    )
    network = Network(
        initial_nodes=['S1_S9', 'S2_T10'],
        resources=resources,
    )

    network.complete_connection(
        maxlen=2,
        algorithm='bfs',
        only_signed=True,
        connect_with_bias=True,
    )

    assert set(network.nodes['Uniprot']) >= {'S1_S9', 'S1', 'S2_T10'}
    assert {
        ('S1_S9', 'S1'),
        ('S1', 'S2_T10'),
    } <= set(network.edges[['source', 'target']].itertuples(
        index=False,
        name=None,
    ))
    assert 'S1_S9' not in calls
    assert 'S2_T10' not in calls
