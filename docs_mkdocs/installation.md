# Installation

## Requirements

NeKo requires **Python ≥ 3.10**.

---

## Install from PyPI (recommended)

```bash
pip install nekomata
```

!!! warning "Package name"
    The PyPI package is called **nekomata** — `pip install neko` and `pip install pyneko` are unrelated packages.

---

## Install from source

For the latest development version:

```bash
git clone https://github.com/sysbio-curie/Neko.git
cd Neko
pip install .
```

---

## Optional: install with Poetry

```bash
git clone https://github.com/sysbio-curie/Neko.git
cd Neko
poetry install
```

---

## System dependencies

### Graphviz

Graphviz is required for network rendering. Install it with your system package manager before installing NeKo.

=== "Linux (apt)"

    ```bash
    sudo apt-get install python3-dev graphviz libgraphviz-dev
    ```

=== "macOS (Homebrew)"

    ```bash
    brew install graphviz
    ```

=== "Windows (Chocolatey)"

    ```bash
    choco install graphviz
    ```

More information: <https://graphviz.org/download/>

### Pandoc (optional, for notebook export)

```bash
sudo apt-get install pandoc   # Linux
brew install pandoc           # macOS
```

---

## Verify installation

```python
import neko
print(neko.__version__)
```

---

## Troubleshooting

### `graphviz` import errors

If `import graphviz` fails, ensure the system Graphviz binary is on your `PATH`:

```bash
which dot   # should point to the graphviz dot binary
```

### OmniPath and UniProt connectivity

NeKo downloads interaction data from OmniPath when building an OmniPath
universe. Identifier translation lazily downloads a reviewed-human mapping
table from UniProt and caches it locally. SIGNOR similarly caches its validated
interaction table and entity dictionaries, so later loads can run offline. If
a resource contains ChEBI accessions, NeKo also caches the compressed ChEBI
compounds table and a small accession-to-name subset. The names are optional:
if ChEBI is unreachable, the canonical accession is used as the display label
and network construction continues. Ensure the first use has network access,
or provide an existing cache with `NEKO_CACHE_DIR`.

```python
from neko.inputs import Universe

resources = Universe("omnipath")
```
