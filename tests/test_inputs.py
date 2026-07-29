import os

import pytest

import pandas as pd

from neko.inputs._db.omnipath import omnipath_universe
import neko.inputs._db.signor as signor_module
import neko.inputs._universe as universe_module

SIGNOR_TSV = (
    b'IDA\tIDB\tDIRECT\tEFFECT\tANNOTATOR\tPMID\tSIGNOR_ID\n'
    b'SRC\tTP53\tt\tup-regulates activity\tlperfetto\t1\tSIGNOR-1\n'
    b'SRC\tTP53\tf\tup-regulates activity\t\t2\tSIGNOR-2\n'
)
INVALID_SIGNOR_RESPONSE = (
    b'<br />\n<b>Fatal error</b>: Allowed memory size exhausted'
)
SIGNOR_ENTITY_CSV = {
    'Download complex data': (
        b'"SIGNOR ID";"COMPLEX NAME";"LIST OF ENTITIES"\n'
        b'SIGNOR-C1;Test complex;"P23511, Q13952"\n'
        b'SIGNOR-C2;Nested complex;"SIGNOR-C1, SIGNOR-PF1"\n'
    ),
    'Download protein family data': (
        b'"SIGNOR ID";"PROT. FAMILY NAME";"LIST OF ENTITIES"\n'
        b'SIGNOR-PF1;ERK1/2;"P27361, P28482"\n'
    ),
    'Download phenotype data': (
        b'"SIGNOR ID";"PHENOTYPE NAME";"PHENOTYPE DESCRIPTION"\n'
        b'SIGNOR-PH1;Cell_death;Programmed cell death\n'
        b'SIGNOR-PH1;Cell_death;Programmed cell death\n'
    ),
    'Download stimulus data': (
        b'"SIGNOR ID";"STIMULUS NAME";"STIMULUS DESCRIPTION"\n'
        b'SIGNOR-ST1;DNA_damage;DNA damage\n'
    ),
}


class FakeResponse:

    def __init__(self, content, error=None):
        self.content = content
        self.error = error

    def raise_for_status(self):
        if self.error:
            raise self.error


@pytest.fixture
def mock_signor_download(monkeypatch, tmp_path):
    calls = {'get': [], 'post': []}

    monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path / 'cache'))

    def get(url, timeout):
        calls['get'].append((url, timeout))
        return FakeResponse(SIGNOR_TSV)

    def post(url, files, timeout):
        submit = files['submit'][1]
        calls['post'].append((url, submit, timeout))
        return FakeResponse(SIGNOR_ENTITY_CSV[submit])

    monkeypatch.setattr(signor_module.requests, 'get', get)
    monkeypatch.setattr(signor_module.requests, 'post', post)

    return calls


def test_download_signor_database(tmp_path, mock_signor_download):
    save_path = tmp_path / 'SIGNOR_Human.tsv'

    signor_module.download_signor_database(str(save_path))

    assert save_path.exists()
    df = pd.read_csv(save_path, sep='\t')
    assert not df.empty

    processed = signor_module.signor(str(save_path))

    assert not processed.empty
    assert {'source', 'target'}.issubset(processed.columns)
    assert processed.loc[0, 'curation_effort'] == 'lperfetto'
    assert bool(processed.loc[0, 'is_directed']) is True
    assert bool(processed.loc[0, 'is_direct']) is True
    assert bool(processed.loc[0, 'consensus_stimulation']) is True
    assert len(mock_signor_download['get']) == 1
    assert len(mock_signor_download['post']) == 4


def test_download_signor_database_as_df(mock_signor_download):
    df = signor_module.download_signor_database()

    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert {'IDA', 'IDB'}.issubset(df.columns)
    assert len(mock_signor_download['get']) == 1
    assert not mock_signor_download['post']


def test_signor_uses_and_populates_managed_database_cache(
    mock_signor_download,
):
    first = signor_module.signor(normalize_entities=False)
    cache_path = signor_module._database_cache_path()

    assert not first.empty
    assert cache_path.exists()
    assert len(mock_signor_download['get']) == 1

    second = signor_module.signor(normalize_entities=False)

    assert first.equals(second)
    assert len(mock_signor_download['get']) == 1


def test_signor_replaces_invalid_managed_database_cache(
    mock_signor_download,
    caplog,
):
    cache_path = signor_module._database_cache_path()
    cache_path.parent.mkdir(parents=True)
    cache_path.write_bytes(
        b'IDA\tIDB\tDIRECT\tEFFECT\tANNOTATOR\tPMID\tSIGNOR_ID\n'
        b'SRC\tTP53\tMAYBE\tup-regulates activity\tcurator\t1\tSIGNOR-1\n'
    )

    processed = signor_module.signor(normalize_entities=False)

    assert not processed.empty
    assert len(mock_signor_download['get']) == 1
    assert 'Ignoring invalid SIGNOR cache' in caplog.text
    cached = pd.read_csv(cache_path, sep='\t')
    assert {'IDA', 'IDB', 'DIRECT'}.issubset(cached.columns)


def test_download_signor_entity_dictionaries(tmp_path, mock_signor_download):
    dictionaries = signor_module.download_signor_entity_dictionaries(tmp_path)

    assert set(dictionaries) == {
        'complexes',
        'protein_families',
        'phenotypes',
        'stimuli',
    }
    assert len(dictionaries['phenotypes']) == 1
    assert len(mock_signor_download['post']) == 4

    for config in signor_module.SIGNOR_ENTITY_DOWNLOADS.values():
        saved = pd.read_csv(tmp_path / config['filename'], sep=';')
        assert not saved.empty


def test_signor_uses_and_populates_managed_entity_cache(
    tmp_path,
    mock_signor_download,
    monkeypatch,
):
    path = tmp_path / 'SIGNOR_Human.tsv'
    path.write_bytes(SIGNOR_TSV)

    first = signor_module.signor(str(path))

    assert not first.empty
    assert len(mock_signor_download['post']) == 4

    monkeypatch.setattr(
        signor_module.requests,
        'post',
        lambda *args, **kwargs: pytest.fail('unexpected dictionary download'),
    )

    second = signor_module.signor(str(path))

    assert first.equals(second)
    assert len(mock_signor_download['post']) == 4

    for config in signor_module.SIGNOR_ENTITY_DOWNLOADS.values():
        assert (signor_module._cache_dir() / config['filename']).exists()


def test_normalize_signor_entities_expands_nested_members():
    dictionaries = {
        entity_type: signor_module._parse_signor_entity_response(
            SIGNOR_ENTITY_CSV[config['submit']],
            config['columns'],
        )
        for entity_type, config in signor_module.SIGNOR_ENTITY_DOWNLOADS.items()
    }
    raw = pd.DataFrame({
        'IDA': [
            'SIGNOR-C1',
            'SIGNOR-C2',
            'SIGNOR-PF1',
            'SIGNOR-PH1',
            'SIGNOR-ST1',
            'SIGNOR-C999',
        ],
        'ENTITYA': [
            'Test complex',
            'Nested complex',
            'ERK1/2',
            'Cell_death',
            'DNA_damage',
            'New complex',
        ],
        'IDB': ['P04637'] * 6,
        'ENTITYB': ['TP53'] * 6,
    })

    normalized = signor_module.normalize_signor_entities(raw, dictionaries)

    assert normalized['IDA'].tolist() == [
        'COMPLEX:P23511_Q13952',
        'COMPLEX:P23511_Q13952_P27361_P28482',
        'PROTEIN_FAMILY:ERK1/2',
        'PHENOTYPE:Cell_death',
        'STIMULUS:DNA_damage',
        'COMPLEX_NAME:New complex',
    ]
    assert raw.loc[0, 'IDA'] == 'SIGNOR-C1'


def test_signor_can_skip_entity_dictionary_download(
    tmp_path,
    monkeypatch,
):
    save_path = tmp_path / 'SIGNOR_Human.tsv'
    save_path.write_bytes(SIGNOR_TSV)
    monkeypatch.setattr(
        signor_module.requests,
        'post',
        lambda *args, **kwargs: pytest.fail('unexpected dictionary download'),
    )

    processed = signor_module.signor(
        str(save_path),
        normalize_entities=False,
    )

    assert processed.loc[0, 'source'] == 'SRC'


def test_public_signor_missing_path_falls_back_with_warning(
    tmp_path,
    monkeypatch,
    caplog,
):
    calls = []
    processed = pd.DataFrame({
        'source': ['SRC'],
        'target': ['TP53'],
        'is_directed': [True],
        'is_stimulation': [True],
        'is_inhibition': [False],
        'form_complex': [False],
    })

    def fake_signor(**kwargs):
        calls.append(kwargs)
        return processed

    monkeypatch.setattr(universe_module._signor, 'signor', fake_signor)
    missing = tmp_path / 'missing.tsv'

    resources = universe_module.signor(
        str(missing),
        normalize_entities=False,
    )

    assert not resources.interactions.empty
    assert calls == [{'normalize_entities': False}]
    assert 'Falling back to the managed NeKo cache' in caplog.text


def test_signor_accepts_provider_directness_values(tmp_path):
    path = tmp_path / 'directness.tsv'
    path.write_bytes(
        b'IDA\tIDB\tDIRECT\tEFFECT\tANNOTATOR\tPMID\tSIGNOR_ID\n'
        b'SRC\tTP53\tt\tup-regulates activity\tcurator\t1\tSIGNOR-1\n'
        b'SRC\tMDM2\tf\tup-regulates activity\tcurator\t2\tSIGNOR-2\n'
    )

    processed = signor_module.signor(
        str(path),
        normalize_entities=False,
    )

    directness = processed.set_index('target')['is_direct'].to_dict()
    assert bool(directness['TP53']) is True
    assert bool(directness['MDM2']) is False


def test_signor_preserves_conflicting_effects_without_false_consensus(
    tmp_path,
):
    path = tmp_path / 'conflicting.tsv'
    path.write_bytes(
        SIGNOR_TSV
        + b'SRC\tTP53\tNO\tdown-regulates activity\tother\t3\tSIGNOR-3\n'
    )

    processed = signor_module.signor(
        str(path),
        normalize_entities=False,
    )

    assert len(processed) == 1
    assert bool(processed.loc[0, 'is_stimulation']) is True
    assert bool(processed.loc[0, 'is_inhibition']) is True
    assert bool(processed.loc[0, 'consensus_stimulation']) is True
    assert bool(processed.loc[0, 'consensus_inhibition']) is False


def test_signor_all_unknown_effects_returns_empty_dataframe(tmp_path):
    path = tmp_path / 'unknown.tsv'
    path.write_bytes(
        b'IDA\tIDB\tDIRECT\tEFFECT\tANNOTATOR\tPMID\tSIGNOR_ID\n'
        b'SRC\tTP53\tYES\tunknown\tcurator\t1\tSIGNOR-1\n'
    )

    processed = signor_module.signor(
        str(path),
        normalize_entities=False,
    )

    assert processed.empty
    assert {'source', 'target', 'consensus_stimulation'}.issubset(
        processed.columns,
    )


def test_signor_rejects_unknown_directness_value(tmp_path):
    path = tmp_path / 'invalid_direct.tsv'
    path.write_bytes(
        b'IDA\tIDB\tDIRECT\tEFFECT\tANNOTATOR\tPMID\tSIGNOR_ID\n'
        b'SRC\tTP53\tMAYBE\tup-regulates activity\tcurator\t1\tSIGNOR-1\n'
    )

    with pytest.raises(ValueError, match='DIRECT column'):
        signor_module.signor(
            str(path),
            normalize_entities=False,
        )


def test_download_signor_database_retries_invalid_response(monkeypatch):
    responses = [
        FakeResponse(INVALID_SIGNOR_RESPONSE),
        FakeResponse(SIGNOR_TSV),
    ]
    delays = []

    monkeypatch.setattr(
        signor_module.requests,
        'get',
        lambda *args, **kwargs: responses.pop(0),
    )
    monkeypatch.setattr(signor_module.time, 'sleep', delays.append)

    df = signor_module.download_signor_database(attempts=2, backoff=0.25)

    assert {'IDA', 'IDB'}.issubset(df.columns)
    assert delays == [0.25]


def test_download_signor_database_rejects_invalid_response(
    tmp_path,
    monkeypatch,
):
    save_path = tmp_path / 'invalid.tsv'
    monkeypatch.setattr(
        signor_module.requests,
        'get',
        lambda *args, **kwargs: FakeResponse(INVALID_SIGNOR_RESPONSE),
    )
    monkeypatch.setattr(signor_module.time, 'sleep', lambda delay: None)

    with pytest.raises(RuntimeError, match='after 2 attempts'):
        signor_module.download_signor_database(
            str(save_path),
            attempts=2,
        )

    assert not save_path.exists()


@pytest.mark.live
@pytest.mark.skipif(
    os.environ.get('NEKO_RUN_LIVE_TESTS') != '1',
    reason='set NEKO_RUN_LIVE_TESTS=1 to call the live SIGNOR service',
)
def test_signor_live(tmp_path, monkeypatch):
    monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path / 'cache'))
    df = signor_module.signor(normalize_entities=False)

    assert {'source', 'target', 'is_direct'}.issubset(df.columns)
    assert not df.empty


def test_omnipath_universe():
    df = omnipath_universe()
    assert isinstance(df, pd.DataFrame)
    assert {'source', 'target'}.issubset(df.columns)
