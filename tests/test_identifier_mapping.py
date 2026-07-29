"""
Unit tests for neko.inputs.identifier_mapping.

These tests are deliberately hermetic (no live network calls): the offline
lookup dicts are seeded directly, and the live-fallback function is
monkeypatched, so behaviour is deterministic and fast. Real end-to-end
translation (including the live download / job-based fallback) is exercised
by the higher-level Network tests in test_network.py / test_strategies.py.
"""

import json

import pytest

from neko.inputs import identifier_mapping as im


@pytest.fixture(autouse=True)
def reset_state(monkeypatch):
    """Ensure each test starts from a clean, isolated module state."""

    monkeypatch.setitem(im._state, 'attempted', True)
    monkeypatch.setitem(im._state, 'symbol_to_uniprot', {})
    monkeypatch.setitem(im._state, 'uniprot_to_symbol', {})
    monkeypatch.setitem(im._state, 'symbol_fallback_cache', {})
    monkeypatch.setitem(im._state, 'uniprot_fallback_cache', {})


def _no_network_fallback(*args, **kwargs):
    raise AssertionError('live fallback should not be called for a cache hit')


class TestLooksLikeUniprotAccession:

    @pytest.mark.parametrize('value', ['P12931', 'Q9Y243', 'A0A0C5B5G6', 'P12931-2'])
    def test_valid_accessions(self, value):
        assert im.looks_like_uniprot_accession(value)

    @pytest.mark.parametrize('value', ['SRC', 'TP53', '', None, 'COMPLEX:A_B'])
    def test_invalid_accessions(self, value):
        assert not im.looks_like_uniprot_accession(value)


class TestToUniprot:

    def test_offline_cache_hit(self, monkeypatch):
        im._state['symbol_to_uniprot']['SRC'] = 'P12931'
        monkeypatch.setattr(im, '_fallback_translate', _no_network_fallback)

        assert im.to_uniprot('SRC') == 'P12931'

    def test_echoes_back_unrecognized_accession(self, monkeypatch):
        monkeypatch.setattr(im, '_fallback_translate', _no_network_fallback)

        assert im.to_uniprot('Q9Y2X3') == 'Q9Y2X3'

    def test_uses_live_fallback_result(self, monkeypatch):
        monkeypatch.setattr(
            im, '_fallback_translate', lambda ids, *a, **k: {'FAKEGENE': 'P00000'},
        )

        assert im.to_uniprot('FAKEGENE') == 'P00000'

    def test_returns_none_when_unmappable(self, monkeypatch):
        monkeypatch.setattr(im, '_fallback_translate', lambda *a, **k: {})

        assert im.to_uniprot('NOT_A_REAL_GENE') is None

    def test_empty_input_returns_none(self):
        assert im.to_uniprot('') is None
        assert im.to_uniprot(None) is None

    def test_fallback_is_memoized_across_calls(self, monkeypatch):
        calls = []

        def _fake_fallback(ids, *a, **k):
            calls.append(set(ids))
            return {}

        monkeypatch.setattr(im, '_fallback_translate', _fake_fallback)

        assert im.to_uniprot('NOT_A_REAL_GENE') is None
        assert im.to_uniprot('NOT_A_REAL_GENE') is None
        assert len(calls) == 1


class TestToGenesymbol:

    def test_offline_cache_hit(self, monkeypatch):
        im._state['uniprot_to_symbol']['P12931'] = 'SRC'
        monkeypatch.setattr(im, '_fallback_translate', _no_network_fallback)

        assert im.to_genesymbol('P12931') == 'SRC'

    def test_echoes_back_when_not_accession_shaped(self, monkeypatch):
        monkeypatch.setattr(im, '_fallback_translate', _no_network_fallback)

        assert im.to_genesymbol('SRC') == 'SRC'

    def test_returns_none_for_unmappable_accession(self, monkeypatch):
        monkeypatch.setattr(im, '_fallback_translate', lambda *a, **k: {})

        assert im.to_genesymbol('Q9Y2X3') is None

    def test_empty_input_returns_none(self):
        assert im.to_genesymbol('') is None
        assert im.to_genesymbol(None) is None


class TestCacheValidation:

    def test_read_cache_missing_files_returns_none(self, tmp_path, monkeypatch):
        monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path))

        assert im._read_cache() is None

    def test_read_cache_rejects_truncated_file(self, tmp_path, monkeypatch):
        monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path))
        cache_dir = im._cache_dir()
        cache_dir.mkdir(parents=True, exist_ok=True)

        table_path = im._table_path()
        table_path.write_text('Entry\tGene Names (primary)\nP12931\tSRC\n')

        with open(im._meta_path(), 'w') as fh:
            json.dump({'rows': 1, 'fetched_at': 0}, fh)

        # Row count (1) is far below _MIN_EXPECTED_ROWS, so the cache must
        # be treated as invalid/untrustworthy rather than silently used.
        assert im._read_cache() is None

    def test_read_cache_accepts_valid_file(self, tmp_path, monkeypatch):
        monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path))
        monkeypatch.setattr(im, '_MIN_EXPECTED_ROWS', 2)
        cache_dir = im._cache_dir()
        cache_dir.mkdir(parents=True, exist_ok=True)

        rows = [('P12931', 'SRC'), ('P04637', 'TP53')]
        im._write_cache(rows)

        assert im._read_cache() == rows

    def test_write_then_read_round_trip(self, tmp_path, monkeypatch):
        monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path))
        monkeypatch.setattr(im, '_MIN_EXPECTED_ROWS', 0)

        rows = [('P12931', 'SRC')]
        im._write_cache(rows)
        loaded = im._read_cache()

        assert loaded == rows

    def test_read_cache_rejects_checksum_mismatch(self, tmp_path, monkeypatch):
        monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path))
        monkeypatch.setattr(im, '_MIN_EXPECTED_ROWS', 1)
        im._write_cache([('P12931', 'SRC')])
        im._table_path().write_text(
            'Entry\tGene Names (primary)\nP04637\tTP53\n',
            encoding='utf-8',
        )

        assert im._read_cache() is None

    def test_read_cache_rejects_non_object_metadata(self, tmp_path, monkeypatch):
        monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path))
        cache_dir = im._cache_dir()
        cache_dir.mkdir(parents=True, exist_ok=True)
        im._table_path().write_text(
            'Entry\tGene Names (primary)\nP12931\tSRC\n',
            encoding='utf-8',
        )
        im._meta_path().write_text('[]', encoding='utf-8')

        assert im._read_cache() is None

    def test_short_download_is_not_cached_or_loaded(self, monkeypatch):
        monkeypatch.setitem(im._state, 'attempted', False)
        monkeypatch.setattr(im, '_read_cache', lambda: None)
        monkeypatch.setattr(
            im,
            '_fetch_reviewed_human_table',
            lambda: [('P12931', 'SRC')],
        )
        writes = []
        monkeypatch.setattr(im, '_write_cache', writes.append)

        im._ensure_loaded()

        assert writes == []
        assert im._state['symbol_to_uniprot'] == {}

    def test_failed_refresh_keeps_loaded_mapping(self, monkeypatch):
        im._state['symbol_to_uniprot']['SRC'] = 'P12931'
        im._state['uniprot_to_symbol']['P12931'] = 'SRC'
        monkeypatch.setattr(im, '_read_cache', lambda: None)
        monkeypatch.setattr(
            im,
            '_fetch_reviewed_human_table',
            lambda: (_ for _ in ()).throw(OSError('offline')),
        )

        im.refresh_cache()

        assert im._state['symbol_to_uniprot'] == {'SRC': 'P12931'}
        assert im._state['uniprot_to_symbol'] == {'P12931': 'SRC'}


class TestBuildDicts:

    def test_build_dicts_basic(self):
        rows = [('P12931', 'SRC'), ('P04637', 'TP53'), ('Q00000', '')]
        symbol_to_uniprot, uniprot_to_symbol = im._build_dicts(rows)

        assert symbol_to_uniprot == {'SRC': 'P12931', 'TP53': 'P04637'}
        assert uniprot_to_symbol == {
            'P12931': 'SRC', 'P04637': 'TP53', 'Q00000': '',
        }

    def test_build_dicts_keeps_first_on_duplicate_symbol(self):
        rows = [('P00001', 'DUP'), ('P00002', 'DUP')]
        symbol_to_uniprot, _ = im._build_dicts(rows)

        assert symbol_to_uniprot['DUP'] == 'P00001'


def test_get_next_link_finds_next_relation_after_other_links():
    headers = {
        'Link': (
            '<https://rest.uniprot.org/first>; rel="first", '
            '<https://rest.uniprot.org/next>; rel="next"'
        ),
    }

    assert im._get_next_link(headers) == 'https://rest.uniprot.org/next'
