# Contributing

We welcome contributions of all kinds — bug reports, documentation fixes, new
features, and new database connectors.

---

## Development setup

```bash
git clone https://github.com/sysbio-curie/Neko.git
cd Neko
poetry install
pre-commit install
```

---

## Code style

NeKo uses the following formatters and linters (enforced via pre-commit):

| Tool | Purpose |
|---|---|
| **black** | Code formatting |
| **autopep8** | PEP 8 compliance |
| **isort** | Import ordering |
| **flake8** | Error / style linting |
| **prettier** | YAML / JSON / Markdown formatting |

Run checks manually:

```bash
pre-commit run --all-files
```

---

## Running tests

```bash
pytest tests/
# or with coverage:
pytest --cov=neko tests/
```

---

## Docstring style

NeKo uses **NumPy-style** docstrings:

```python
def my_function(x: int, y: float) -> str:
    """
    One-line summary.

    Longer description if needed.

    Parameters
    ----------
    x : int
        Description of x.
    y : float
        Description of y.

    Returns
    -------
    str
        Description of the return value.
    """
```

---

## Pull request checklist

- [ ] Tests pass (`pytest tests/`)
- [ ] New code is covered by tests
- [ ] Docstrings follow NumPy style
- [ ] `pre-commit run --all-files` passes
- [ ] Changelog entry added

---

## Reporting issues

Open a [GitHub issue](https://github.com/sysbio-curie/Neko/issues) and include:

1. NeKo version (`python -c "import neko; print(neko.__version__)"`)
2. Python version
3. OS
4. Minimal reproducible example
