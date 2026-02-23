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

### OmniPath / PyPath connectivity

NeKo downloads interaction data from the OmniPath web service on first use.  
Ensure you have a working internet connection when calling `Universe().build()`.
