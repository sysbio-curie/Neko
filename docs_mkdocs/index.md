# NeKo – Network Konstructor

<p align="center">
  <img src="assets/neko_logo.png" width="120" alt="NeKo Logo"/>
</p>

[![Tests](https://github.com/sysbio-curie/Neko/actions/workflows/build.yaml/badge.svg)](https://github.com/sysbio-curie/Neko/actions/workflows/build.yaml)
[![Docs](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://sysbio-curie.github.io/Neko/)
[![PyPI](https://img.shields.io/pypi/v/nekomata.svg)](https://pypi.org/project/nekomata/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**NeKo** is a Python package for extracting, visualising, converting, and studying interactions from databases into executable activity flow-based models. It is built on top of [OmniPath](https://github.com/saezlab/omnipath), [PyPath](https://github.com/saezlab/pypath), and [Atopo](https://github.com/druglogics/atopo).

---

## Key features

| Feature | Description |
|---|---|
| **Network creation** | Build signalling networks from curated prior-knowledge databases |
| **Node connection** | Connect nodes and subnetworks with flexible stepwise strategies |
| **Gene-to-phenotype** | Map gene sets to phenotypic categories via Gene Ontology |
| **Visualisation** | Render networks with Graphviz or the interactive yFiles widget |
| **Export** | Export to SIF, GML, GraphML, BND/CFG logical model formats |
| **Network history** | Automatic snapshots, branching state management, and HTML/SVG rendering |
| **Comparison** | Compare two networks and highlight differences |

---

## Citation

If you use NeKo in your research, please cite:

> Ruscone M, Tsirvouli E, Checcoli A, Turei D, Barillot E, et al. (2025)
> **NeKo: A tool for automatic network construction from prior knowledge.**
> *PLOS Computational Biology* 21(9): e1013300.
> <https://doi.org/10.1371/journal.pcbi.1013300>

---

## Quick start

```python
from neko.core.network import Network
from neko.inputs import Universe
from neko._visual.visualize_network import NetworkVisualizer

# Load the interaction universe from OmniPath
resources = Universe()
resources.build()

# Create a network from a list of genes
net = Network(["EGFR", "KRAS", "TP53", "AKT1"], resources=resources.interactions)

# Connect nodes using the available database interactions
net.connect_nodes()

# Visualise
vis = NetworkVisualizer(net)
vis.render()
```

See the [Tutorials](tutorials/index.md) for full worked examples, or jump to the [API Reference](api/index.md) for detailed documentation.

---

## Comparison with Sphinx docs

During the transition period both documentation backends are live. Feel free to use whichever you prefer and let us know via a [GitHub issue](https://github.com/sysbio-curie/Neko/issues):

- **This site (MkDocs / Material)** – <https://sysbio-curie.github.io/Neko/>
- **Sphinx / RTD** – <https://sysbio-curie.github.io/Neko/sphinx/>
