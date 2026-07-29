import pandas as pd
import pytest

from neko.core.network import Network
from neko.core import tools
from neko.inputs import chebi_mapping


@pytest.fixture(autouse=True)
def isolated_chebi_cache(tmp_path, monkeypatch):
    monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path / 'cache'))
    monkeypatch.setitem(chebi_mapping._state, 'loaded', False)
    monkeypatch.setitem(chebi_mapping._state, 'id_to_name', {})
    monkeypatch.setitem(chebi_mapping._state, 'missing', set())


def test_normalize_identifier():
    assert chebi_mapping.normalize_identifier('CHEBI:16618') == 'CHEBI:16618'
    assert chebi_mapping.normalize_identifier('chebi_16618') == 'CHEBI:16618'
    assert chebi_mapping.normalize_identifier('P12931') is None
    assert chebi_mapping.normalize_identifier(None) is None


def test_extract_names_streams_only_requested_accessions(tmp_path, monkeypatch):
    compounds_path = tmp_path / 'compounds.tsv.gz'
    pd.DataFrame({
        'name': ['HTML name', 'Fallback name', 'Unused'],
        'chebi_accession': ['CHEBI:16618', 'CHEBI:1', 'CHEBI:2'],
        'ascii_name': ['ASCII name', pd.NA, 'Unused ASCII'],
    }).to_csv(
        compounds_path,
        sep='\t',
        index=False,
        compression='gzip',
    )
    monkeypatch.setattr(chebi_mapping, '_MIN_EXPECTED_ROWS', 1)

    names = chebi_mapping._extract_names(
        compounds_path,
        {'CHEBI:16618', 'CHEBI:1'},
    )

    assert names == {
        'CHEBI:16618': 'ASCII name',
        'CHEBI:1': 'Fallback name',
    }


def test_ensure_names_uses_derived_cache_without_downloading(monkeypatch):
    chebi_mapping._write_mapping_cache({'CHEBI:16618': 'Cached name'})

    def unexpected_download(*args, **kwargs):
        pytest.fail('a derived-cache hit must not download ChEBI data')

    monkeypatch.setattr(
        chebi_mapping,
        '_download_and_extract',
        unexpected_download,
    )

    assert chebi_mapping.ensure_names({'CHEBI:16618'}) == {
        'CHEBI:16618': 'Cached name',
    }


def test_failed_name_download_falls_back_to_accession(monkeypatch):
    def failed_download(*args, **kwargs):
        raise OSError('offline')

    monkeypatch.setattr(
        chebi_mapping,
        '_download_and_extract',
        failed_download,
    )
    monkeypatch.setattr(
        tools.mapping,
        'to_uniprot',
        lambda value: pytest.fail('ChEBI must not be sent to UniProt'),
    )

    assert tools.mapping_node_identifier('CHEBI:16618') == [
        None,
        'CHEBI:16618',
        'CHEBI:16618',
    ]


@pytest.mark.parametrize(
    'identifier',
    ['CID:53396311', 'URS00005C2A6D_9606', 'SIGNOR-FP6'],
)
def test_other_signor_nonprotein_ids_bypass_uniprot(identifier, monkeypatch):
    monkeypatch.setattr(
        tools.mapping,
        'to_uniprot',
        lambda value: pytest.fail('non-protein IDs must not reach UniProt'),
    )

    assert tools.mapping_node_identifier(identifier) == [
        None,
        identifier,
        identifier,
    ]


def test_connectivity_uses_label_for_legacy_null_identifier():
    from neko.core.tools import is_connected

    network = type('LegacyNetwork', (), {})()
    network.nodes = pd.DataFrame({
        'Genesymbol': ['P1', 'CHEBI:16618'],
        'Uniprot': ['P1', None],
    })
    network.edges = pd.DataFrame({
        'source': ['CHEBI:16618'],
        'target': ['P1'],
    })

    assert is_connected(network)


def test_connect_as_atopo_with_chebi_neighbour_terminates(
    monkeypatch,
    caplog,
):
    names = {'CHEBI:16618': 'Phosphatidylinositol trisphosphate'}
    monkeypatch.setattr(
        chebi_mapping,
        'ensure_names',
        lambda identifiers: {
            identifier: names[identifier]
            for identifier in identifiers
            if identifier in names
        },
    )
    monkeypatch.setattr(chebi_mapping, 'to_name', names.get)
    def protein_mapping(value):
        return None if value == 'UNKNOWN' else value

    monkeypatch.setattr(tools.mapping, 'to_uniprot', protein_mapping)
    monkeypatch.setattr(tools.mapping, 'to_genesymbol', protein_mapping)
    resources = pd.DataFrame({
        'source': ['CHEBI:16618', 'P1'],
        'target': ['P1', 'P2'],
        'is_directed': [True, True],
        'is_stimulation': [True, True],
        'is_inhibition': [False, False],
        'form_complex': [False, False],
    })
    network = Network(initial_nodes=['P1'], resources=resources)

    network.connect_as_atopo(
        strategy='radial',
        max_len=1,
        outputs=['P2', 'UNKNOWN'],
    )

    assert set(network.nodes['Uniprot']) == {'P1', 'P2'}
    assert set(network.edges[['source', 'target']].itertuples(
        index=False,
        name=None,
    )) == {('P1', 'P2')}
    assert 'Ignoring output nodes' in caplog.text


def test_remove_node_accepts_chebi_display_name(monkeypatch):
    label = 'Phosphatidylinositol trisphosphate'
    monkeypatch.setattr(
        chebi_mapping,
        'ensure_names',
        lambda identifiers: {'CHEBI:16618': label},
    )
    monkeypatch.setattr(chebi_mapping, 'to_name', lambda identifier: label)
    monkeypatch.setattr(
        tools.mapping,
        'to_uniprot',
        lambda value: (
            value
            if value == 'P1'
            else pytest.fail('a known display label must not be mapped')
        ),
    )
    monkeypatch.setattr(tools.mapping, 'to_genesymbol', lambda value: value)
    resources = pd.DataFrame({
        'source': ['CHEBI:16618'],
        'target': ['P1'],
        'is_directed': [True],
        'is_stimulation': [True],
        'is_inhibition': [False],
        'form_complex': [False],
    })
    network = Network(initial_nodes=['P1'], resources=resources)
    network.add_node('CHEBI:16618')
    network.add_edge(resources)
    network._history_enabled = False

    network.remove_node(label)

    assert 'CHEBI:16618' not in set(network.nodes['Uniprot'])
    assert network.edges.empty
