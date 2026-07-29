"""Shared helpers for NeKo's on-disk caches."""

from pathlib import Path
import os


def cache_dir(subdirectory: str) -> Path:
    """Return a resource-specific directory below NeKo's cache root."""

    override = os.environ.get('NEKO_CACHE_DIR')

    if override:
        base = Path(override)
    elif os.name == 'nt':
        base = Path(
            os.environ.get(
                'LOCALAPPDATA',
                Path.home() / 'AppData' / 'Local',
            ),
        )
    else:
        base = Path(
            os.environ.get('XDG_CACHE_HOME', Path.home() / '.cache'),
        )

    return base / 'neko' / subdirectory
