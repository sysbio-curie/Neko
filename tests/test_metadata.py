from pathlib import Path

import toml

from neko import _metadata


def test_get_metadata_uses_source_tree_pyproject():
    project = toml.load(Path(__file__).parents[1] / 'pyproject.toml')

    metadata = _metadata.get_metadata()

    assert metadata['name'] == project['tool']['poetry']['name']
    assert metadata['version'] == project['tool']['poetry']['version']


def test_get_metadata_uses_nekomata_distribution(monkeypatch):
    requested_distributions = []

    monkeypatch.setattr(_metadata.os.path, 'exists', lambda path: False)

    def distribution_metadata(name):
        requested_distributions.append(name)
        return {
            'Name': 'nekomata',
            'Version': '1.9.0',
            'License': 'GPL-3.0-only',
        }

    monkeypatch.setattr(
        _metadata.importlib.metadata,
        'metadata',
        distribution_metadata,
    )

    metadata = _metadata.get_metadata()

    assert requested_distributions == ['nekomata']
    assert metadata['name'] == 'nekomata'
    assert metadata['version'] == '1.9.0'


def test_get_metadata_has_explicit_unknown_fallback(monkeypatch):
    monkeypatch.setattr(_metadata.os.path, 'exists', lambda path: False)

    def missing_distribution(name):
        raise _metadata.importlib.metadata.PackageNotFoundError(name)

    monkeypatch.setattr(
        _metadata.importlib.metadata,
        'metadata',
        missing_distribution,
    )

    metadata = _metadata.get_metadata()

    assert metadata['version'] == '0+unknown'
