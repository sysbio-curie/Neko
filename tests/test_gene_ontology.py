import pandas as pd
import pytest
import requests

from neko._annotations import gene_ontology
from neko._annotations.gene_ontology import (
    AnnotationServiceError,
    GeneOntologyError,
    GeneOntologyNotFoundError,
    Ontology,
)


class _Response:
    def __init__(self, payload=None, status_code=200, json_error=None):
        self.payload = payload
        self.status_code = status_code
        self.json_error = json_error

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def json(self):
        if self.json_error:
            raise self.json_error
        return self.payload


class _Session:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = []

    def get(self, url, params=None, timeout=None):
        self.calls.append({
            'url': url,
            'params': params,
            'timeout': timeout,
        })
        return self.responses.pop(0)


def _term():
    return {
        'goid': 'GO:0062043',
        'label': (
            'positive regulation of cardiac epithelial to mesenchymal '
            'transition'
        ),
    }


def _association(
        association_id,
        gene_id,
        symbol,
        *,
        object_id='GO:0062043',
        taxon_id='NCBITaxon:9606',
        negated=False,
        evidence='ECO:0000250',
    ):
    return {
        'id': association_id,
        'subject': {
            'id': gene_id,
            'label': symbol,
            'taxon': {
                'id': taxon_id,
                'label': 'Homo sapiens',
            },
        },
        'object': {'id': object_id},
        'negated': negated,
        'evidence': evidence,
        'evidence_types': [{'id': evidence}],
    }


def test_go_api_paginates_without_num_found_and_deduplicates():
    session = _Session([
        _Response(_term()),
        _Response({'associations': [
            _association('a1', 'UniProtKB:Q04771', 'ACVR1'),
            _association('a2', 'UniProtKB:Q04771', 'ACVR1'),
        ]}),
        _Response({'associations': [
            _association('a3', 'UniProtKB:P46531', 'NOTCH1'),
        ]}),
    ])
    ontology = Ontology(session=session, timeout=12)

    genes = ontology.fetch_go_genes('0062043', page_size=2)

    assert [(gene.gene_id, gene.symbol) for gene in genes] == [
        ('UniProtKB:Q04771', 'ACVR1'),
        ('UniProtKB:P46531', 'NOTCH1'),
    ]
    association_calls = session.calls[1:]
    assert [call['params']['start'] for call in association_calls] == [0, 2]
    assert all(
        call['params']['taxon'] == 'NCBITaxon:9606'
        for call in association_calls
    )
    assert all(call['timeout'] == 12 for call in session.calls)


def test_go_api_applies_filters_locally():
    associations = [
        _association('exact', 'UniProtKB:Q04771', 'ACVR1'),
        _association(
            'descendant',
            'UniProtKB:P17813',
            'ENG',
            object_id='GO:1905007',
        ),
        _association(
            'automatic',
            'UniProtKB:P61812',
            'TGFB2',
            evidence='ECO:0000501',
        ),
        _association(
            'negated',
            'UniProtKB:P46531',
            'NOTCH1',
            negated=True,
        ),
        _association(
            'mouse',
            'UniProtKB:P12345',
            'MouseGene',
            taxon_id='NCBITaxon:10090',
        ),
    ]
    session = _Session([
        _Response(_term()),
        _Response({'associations': associations}),
    ])
    ontology = Ontology(session=session)

    genes = ontology.fetch_go_genes(
        'GO:0062043',
        include_descendants=True,
        exclude_automatic_assertions=True,
    )

    assert [gene.symbol for gene in genes] == ['ACVR1', 'ENG']


def test_go_api_defaults_to_exact_term_associations():
    session = _Session([
        _Response(_term()),
        _Response({'associations': [
            _association('exact', 'UniProtKB:Q04771', 'ACVR1'),
            _association(
                'descendant',
                'UniProtKB:P17813',
                'ENG',
                object_id='GO:1905007',
            ),
        ]}),
    ])

    markers = Ontology(session=session).get_markers(
        id_accession='GO:0062043',
    )

    assert markers == ['ACVR1']


def test_go_api_raises_for_unknown_term():
    ontology = Ontology(session=_Session([_Response(status_code=404)]))

    with pytest.raises(GeneOntologyNotFoundError, match='GO:0062043'):
        ontology.get_term('GO:0062043')


def test_go_api_wraps_invalid_json():
    ontology = Ontology(session=_Session([
        _Response(json_error=ValueError('invalid JSON')),
    ]))

    with pytest.raises(GeneOntologyError, match='invalid JSON'):
        ontology.get_term('GO:0062043')


def test_go_api_wraps_http_errors():
    ontology = Ontology(session=_Session([_Response(status_code=503)]))

    with pytest.raises(GeneOntologyError, match='request failed') as error:
        ontology.get_term('GO:0062043')

    assert isinstance(error.value.__cause__, requests.HTTPError)


def test_go_api_rejects_invalid_association_schema():
    session = _Session([
        _Response(_term()),
        _Response({'associations': None}),
    ])

    with pytest.raises(GeneOntologyError, match="'associations'.*list"):
        Ontology(session=session).fetch_go_genes('GO:0062043')


@pytest.mark.parametrize(
    'go_id',
    [None, '', 'GO:123', 'HP:0062043', 'GO:ABCDEFG'],
)
def test_go_api_validates_accessions(go_id):
    with pytest.raises(ValueError, match='GO accession'):
        Ontology(session=_Session([])).get_term(go_id)


def _annotations():
    return pd.DataFrame([
        {
            'genesymbol': 'SRC',
            'record_id': 1,
            'label': 'tissue',
            'value': ' Colon ',
        },
        {
            'genesymbol': 'SRC',
            'record_id': 1,
            'label': 'level',
            'value': 'Medium',
        },
        {
            'genesymbol': 'TP53',
            'record_id': 2,
            'label': 'tissue',
            'value': 'colon',
        },
        {
            'genesymbol': 'TP53',
            'record_id': 2,
            'label': 'level',
            'value': 'Not detected',
        },
        {
            'genesymbol': 'TP53',
            'record_id': 3,
            'label': 'tissue',
            'value': 'liver cancer',
        },
        {
            'genesymbol': 'TP53',
            'record_id': 3,
            'label': 'level',
            'value': 'High',
        },
    ])


def test_non_cancer_annotations_use_one_omnipath_request(monkeypatch):
    calls = []

    def get(**kwargs):
        calls.append(kwargs)
        return _annotations()

    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        get,
    )
    genes = pd.DataFrame({
        'Genesymbol': ['SRC', 'TP53', 'UNKNOWN', 'SRC'],
    })

    result = Ontology().check_tissue_annotations(
        genes,
        '  COLON ',
    )

    assert calls == [{
        'proteins': ['SRC', 'TP53', 'UNKNOWN'],
        'resources': 'HPA_tissue',
    }]
    assert result.to_dict('records') == [
        {'Genesymbol': 'SRC', 'in_tissue': True},
        {'Genesymbol': 'TP53', 'in_tissue': False},
        {'Genesymbol': 'UNKNOWN', 'in_tissue': False},
        {'Genesymbol': 'SRC', 'in_tissue': True},
    ]


def test_tissue_match_is_not_a_substring(monkeypatch):
    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        lambda **_: _annotations(),
    )

    result = Ontology().check_tissue_annotations(
        pd.DataFrame({'Genesymbol': ['SRC']}),
        'col',
    )

    assert result['in_tissue'].tolist() == [False]


def test_cancer_annotations_use_direct_hpa_data(monkeypatch):
    table = pd.DataFrame([
        {
            'Gene name': 'SRC',
            'Cancer': 'colorectal cancer',
            'High': 0,
            'Medium': 2,
            'Low': 0,
            'Not detected': 1,
        },
        {
            'Gene name': 'TP53',
            'Cancer': 'colorectal cancer',
            'High': 0,
            'Medium': 0,
            'Low': 0,
            'Not detected': 3,
        },
    ])
    monkeypatch.setattr(
        gene_ontology,
        '_load_hpa_cancer_table',
        lambda: table,
    )
    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        lambda **_: pytest.fail('Cancer data must bypass OmniPath'),
    )

    result = Ontology().check_tissue_annotations(
        pd.DataFrame({'Genesymbol': ['SRC', 'TP53', 'UNKNOWN']}),
        ' Colorectal  Cancer ',
    )

    assert result.to_dict('records') == [
        {'Genesymbol': 'SRC', 'in_tissue': True},
        {'Genesymbol': 'TP53', 'in_tissue': False},
        {'Genesymbol': 'UNKNOWN', 'in_tissue': False},
    ]


def test_hpa_cancer_loader_reuses_valid_cache(tmp_path, monkeypatch):
    monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path / 'cache'))
    monkeypatch.setattr(gene_ontology, '_MIN_HPA_CANCER_ROWS', 1)
    path = gene_ontology._hpa_cancer_cache_path()
    path.parent.mkdir(parents=True)
    pd.DataFrame([{
        'Gene name': 'SRC',
        'Cancer': 'colorectal cancer',
        'High': 0,
        'Medium': 1,
        'Low': 0,
        'Not detected': 0,
    }]).to_csv(
        path,
        sep='\t',
        index=False,
        compression={
            'method': 'zip',
            'archive_name': 'cancer_data.tsv',
        },
    )
    monkeypatch.setattr(
        gene_ontology,
        '_download_hpa_cancer_table',
        lambda: pytest.fail('A valid cache must not be downloaded again'),
    )

    table = gene_ontology._load_hpa_cancer_table()

    assert table['Gene name'].tolist() == ['SRC']


def test_tissue_annotations_wraps_download_failure(monkeypatch):
    def fail(**_):
        raise RuntimeError('No active exception to reraise')

    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        fail,
    )

    with pytest.raises(AnnotationServiceError, match='temporarily unavailable') as error:
        Ontology().check_tissue_annotations(
            pd.DataFrame({'Genesymbol': ['SRC']}),
            'colon',
        )

    assert isinstance(error.value.__cause__, RuntimeError)


def test_tissue_annotations_rejects_invalid_response(monkeypatch):
    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        lambda **_: pd.DataFrame({'genesymbol': ['SRC']}),
    )

    with pytest.raises(AnnotationServiceError, match='missing columns'):
        Ontology().check_tissue_annotations(
            pd.DataFrame({'Genesymbol': ['SRC']}),
            'colon',
        )


def test_tissue_annotations_accepts_empty_gene_dataframe(monkeypatch):
    def fail(**_):
        pytest.fail('An empty input must not call OmniPath')

    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        fail,
    )

    result = Ontology().check_tissue_annotations(
        pd.DataFrame({'Genesymbol': []}),
        'colorectal cancer',
    )

    assert result.empty
    assert result.dtypes.to_dict() == {
        'Genesymbol': object,
        'in_tissue': bool,
    }


@pytest.mark.parametrize(
    ('genes', 'error'),
    [
        (pd.DataFrame(), "must contain a 'Genesymbol' column"),
        (pd.DataFrame({'Genesymbol': [None]}), 'missing gene symbol'),
        (pd.DataFrame({'Genesymbol': ['  ']}), 'empty gene symbol'),
    ],
)
def test_tissue_annotations_validates_gene_symbols(genes, error):
    with pytest.raises(ValueError, match=error):
        Ontology().check_tissue_annotations(genes, 'colorectal cancer')


@pytest.mark.parametrize('tissue', [None, '', '  '])
def test_tissue_annotations_validates_tissue(tissue):
    with pytest.raises(ValueError, match='non-empty string'):
        Ontology().check_tissue_annotations(
            pd.DataFrame({'Genesymbol': ['SRC']}),
            tissue,
        )
