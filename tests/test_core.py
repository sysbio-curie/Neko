import pandas as pd
from pypath.utils import mapping


def identity(x):
    return x


def setup_module(module):
    module.orig_id = mapping.id_from_label0
    module.orig_label = mapping.label
    mapping.id_from_label0 = identity
    mapping.label = identity


def teardown_module(module):
    mapping.id_from_label0 = module.orig_id
    mapping.label = module.orig_label


def sample_df():
    return pd.DataFrame({
        'source': ['P1', 'P2', 'P3'],
        'target': ['P2', 'P3', 'P1'],
        'is_directed': [True, True, True],
        'is_stimulation': [True, False, True],
        'is_inhibition': [False, True, False],
        'form_complex': [False, False, False],
    })


def test_universe_basic():
    from neko.inputs._universe import Universe
    df = sample_df()
    u = Universe(df)
    assert len(u) == 3
    assert {'P1', 'P2', 'P3'} <= u.nodes


def test_network_add_edge():
    from neko.core.network import Network
    df = sample_df()
    net = Network(initial_nodes=['P1'], resources=df)
    edge = pd.DataFrame({
        'source': ['P1'],
        'target': ['P2'],
        'type': ['activation'],
        'references': ['ref'],
        'is_stimulation': [True],
        'is_inhibition': [False],
    })
    net.add_edge(edge)
    assert len(net.edges) == 1
    assert {'P1', 'P2'} <= set(net.edges[['source', 'target']].stack())


def test_exports(tmp_path):
    from neko.core.network import Network
    from neko._outputs.exports import Exports
    df = sample_df()
    net = Network(initial_nodes=['P1'], resources=df)
    edge = pd.DataFrame({
        'source': ['P1'],
        'target': ['P2'],
        'type': ['activation'],
        'references': ['ref'],
        'is_stimulation': [True],
        'is_inhibition': [False],
    })
    net.add_edge(edge)
    exp = Exports(net)
    sif_file = tmp_path / 'test.sif'
    bnet_file = tmp_path / 'test.bnet'
    exp.export_sif(str(sif_file))
    exp.export_bnet(str(bnet_file))
    assert sif_file.exists()
    created = list(tmp_path.glob('test*.bnet'))
    assert created

