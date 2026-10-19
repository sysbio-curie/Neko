"""Offline regressions for SIF labels, canonical IDs, and Boolean export."""

import pandas as pd
import pytest

from neko.core.network import Network
from neko._outputs.exports import Exports
from neko.inputs import identifier_mapping as im


@pytest.fixture(autouse=True)
def offline_identifiers(monkeypatch):
    # Exercise the real complex/label mapper with a small local protein table.
    monkeypatch.setitem(im._state, 'attempted', True)
    monkeypatch.setitem(im._state, 'symbol_to_uniprot', {
        'SRC': 'P12931', 'TP53': 'P04637',
    })
    monkeypatch.setitem(im._state, 'uniprot_to_symbol', {
        'P12931': 'SRC', 'P04637': 'TP53',
    })
    monkeypatch.setitem(im._state, 'symbol_fallback_cache', {})
    monkeypatch.setitem(im._state, 'uniprot_fallback_cache', {})
    monkeypatch.setattr(im, '_fallback_translate', lambda *args, **kwargs: {})


def empty_resources():
    return pd.DataFrame(columns=[
        'source', 'target', 'is_directed', 'is_stimulation',
        'is_inhibition', 'form_complex',
    ])

@pytest.mark.parametrize('initial_nodes', [None, []])
@pytest.mark.parametrize(('identifier', 'label'), [
    ('COMPLEX:P12931_P04637', 'COMPLEX:SRC_TP53'),
    ('COMPLEX:unknown_member', 'COMPLEX:unknown_member'),
    ('custom-label', 'custom-label'),
    ('Q9Y2X3', 'Q9Y2X3'),  # Valid accession without a known symbol.
    ('PROTEIN_FAMILY:ERK1/2', 'PROTEIN_FAMILY:ERK1/2'),
    ('P12931', 'SRC'),
])
def test_sif_node_representations_agree(identifier, label, initial_nodes):
    net = Network(initial_nodes=initial_nodes, resources=empty_resources())
    assert net.add_node(identifier, from_sif=True)
    assert net.add_node(identifier, from_sif=True)

    assert net.nodes.to_dict('records') == [{
        'Genesymbol': label, 'Uniprot': identifier, 'Type': 'NaN',
    }]
    assert net.initial_nodes == [label]
    node, = net.nodes_as_objects()
    assert (node.id, node.metadata['Genesymbol'], node.type) == (
        identifier, label, 'NaN',
    )
    before = net.nodes.copy(deep=True)
    net.sync_nodes_to_df()
    pd.testing.assert_frame_equal(net.nodes, before, check_dtype=False)


def test_unmapped_sif_node_is_initial_node(monkeypatch):
    import neko.core.network as network_module

    net = Network(initial_nodes=[], resources=empty_resources())
    monkeypatch.setattr(
        network_module, 'mapping_node_identifier', lambda node: [None] * 3,
    )
    net.add_node('unmapped', from_sif=True)
    assert net.initial_nodes == ['unmapped']
    assert net.nodes.iloc[0]['Genesymbol'] == 'unmapped'
    assert net.nodes.iloc[0]['Uniprot'] == 'unmapped'
    node, = net.nodes_as_objects()
    assert (node.id, node.metadata['Genesymbol']) == ('unmapped', 'unmapped')


def test_sif_to_bnet_preserves_complex_logic_and_network(tmp_path):
    sif = tmp_path / 'complex.sif'
    sif.write_text(
        'P12931 form_complex COMPLEX:P12931_P04637 assembly\n'
        'P04637 form_complex COMPLEX:P12931_P04637 assembly\n'
        'COMPLEX:P12931_P04637 stimulation custom-label regulation\n'
        'P04637 inhibition custom-label regulation\n'
    )
    net = Network(sif_file=str(sif), resources=empty_resources())
    nodes_before = net.nodes.copy(deep=True)
    edges_before = net.edges.copy(deep=True)
    initial_before = net.initial_nodes.copy()
    expected_edges = {
        ('P12931', 'COMPLEX:P12931_P04637', 'assembly', 'form complex', 'SIF file'),
        ('P04637', 'COMPLEX:P12931_P04637', 'assembly', 'form complex', 'SIF file'),
        ('COMPLEX:P12931_P04637', 'custom-label', 'regulation', 'stimulation', 'SIF file'),
        ('P04637', 'custom-label', 'regulation', 'inhibition', 'SIF file'),
    }
    assert set(net.edges.itertuples(index=False, name=None)) == expected_edges
    assert {
        (edge.source, edge.target, edge.interaction_type,
         edge.metadata['Effect'], edge.evidence)
        for edge in net.edges_as_objects()
    } == expected_edges
    assert set(net.initial_nodes) == {'SRC', 'TP53', 'COMPLEX:SRC_TP53', 'custom-label'}
    assert {node.id: node.metadata['Genesymbol'] for node in net.nodes_as_objects()} == {
        'P12931': 'SRC', 'P04637': 'TP53',
        'COMPLEX:P12931_P04637': 'COMPLEX:SRC_TP53',
        'custom-label': 'custom-label',
    }

    exporter = Exports(net)
    assert exporter.interactions['References'].tolist() == ['SIF file'] * 4
    exporter.export_bnet(str(tmp_path / 'complex.bnet'))
    rules = dict(line.split(', ', 1) for line in
                 (tmp_path / 'complex_1.bnet').read_text().splitlines()[2:])
    assert rules == {
        'SRC': 'SRC', 'TP53': 'TP53',
        'COMPLEX_SRC_TP53': '(SRC & TP53)',
        'custom_label': '(COMPLEX_SRC_TP53) & !(TP53)',
    }
    pd.testing.assert_frame_equal(net.nodes, nodes_before)
    pd.testing.assert_frame_equal(net.edges, edges_before)
    assert net.initial_nodes == initial_before


def test_sif_complex_sanitization_collision_rejected(tmp_path):
    sif = tmp_path / 'collision.sif'
    sif.write_text('COMPLEX:P12931_P04637 stimulation COMPLEX_SRC_TP53\n')
    net = Network(sif_file=str(sif), resources=empty_resources())
    with pytest.raises(ValueError, match='collide after BNet identifier sanitization'):
        Exports(net).export_bnet(str(tmp_path / 'collision.bnet'))
    assert not list(tmp_path.glob('*.bnet'))
