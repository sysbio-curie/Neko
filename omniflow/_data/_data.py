from typing import Any, Callable
import os
import pathlib as pl

import pypath_common.misc as _common

"""
Built-in data of the module.

To access the data,
"""

_DATADIR = pl.Path(__file__).parent.absolute()
_FORMATS = {
    'json': 'json.load',
    'yaml': 'yaml.safe_load',
    'csv': 'pandas.read_csv',
    'tsv': 'pandas.read_csv',
}


def load(label: str, reader: Callable | None, **kwargs) -> Any:
    """
    Load data shipped with the module or data from a path.

    Args:
        label:
            Label of a built-in dataset or path to a file.

    Returns:
        The object read from the file (typically a dict or list).
    """

    if os.path.exists(label):

        path = pl.Path(label)

    else:

        path = None

        for d, dirs, files in os.walk(_DATADIR):

            for f in files:

                p = pl.Path(d) / f

                if (
                    label == p.stem or
                    label == p.name or
                    label == str(p).rsplit('.', maxsplit = 1)[0]
                ):

                    path = p
                    break

            if path:

                break

    if path and os.path.exists(path):

        path = pl.Path(path).absolute()

        if not reader:

            ext = path.name.rsplit('.', maxsplit = 1)[-1].lower()
            if ext == 'tsv':
                kwargs['sep'] = '\t'
            reader = _FORMATS.get(ext, lambda x: x.readlines())

        if not callable(reader):

            reader = _misc.from_module(reader)

        with open(path) as fp:

            return reader(fp, **kwargs)


def builtins() -> list[str]:
    """
    List of built-in datasets.
    """

    return [
        f.rsplit('.', maxsplit = 1)[0]
        for d, dirs, files in os.walk(_DATADIR):
        for f in files
        if f.rsplit('.', maxsplit = 1)[-1].lower() in _FORMATS
    ]


