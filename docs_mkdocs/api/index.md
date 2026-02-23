# API Reference

NeKo's public API is organised into five main areas:

| Module | Description |
|---|---|
| [`neko.core.network.Network`](network.md) | Core class – build, manipulate, and query a signalling network |
| [`neko._visual.visualize_network.NetworkVisualizer`](visualizer.md) | Render networks with Graphviz or the yFiles widget |
| [`neko._methods.enrichment_methods.Connections`](connections.md) | Algorithms for enriching a network from an interaction database |
| [`neko._annotations.gene_ontology.Ontology`](ontology.md) | Gene Ontology utilities and phenotype mapping |
| [`neko._outputs.exports`](exports.md) | Export helpers (SIF, GML, GraphML, BND/CFG) |

---

## Design philosophy

NeKo follows a **Network-centric** design:

1. Start with a list of gene/protein identifiers.
2. Attach an interaction `Universe` from OmniPath (or custom CSV/DataFrame).
3. Use `Network` methods to connect, expand, and annotate nodes.
4. Visualise or export the result.

All mutation methods on `Network` automatically create a snapshot in the branching **NetworkHistory**, so every intermediate state is recoverable.

---

## Import conventions

```python
# Core
from neko.core.network import Network

# Visualisation
from neko._visual.visualize_network import NetworkVisualizer

# Interaction universe
from neko.inputs import Universe

# Ontology
from neko._annotations.gene_ontology import Ontology
```
