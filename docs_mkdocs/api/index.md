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

## SIGNOR entities

`neko.inputs.signor()` normalizes SIGNOR-specific nodes automatically. It uses
NeKo's validated local cache for the interaction table and the complex,
protein-family, phenotype, and stimulus dictionaries, downloading only the
missing resources. It expands complexes recursively into OmniPath-compatible
`COMPLEX:` identifiers and assigns readable typed identifiers to the other
entity classes.

```python
from neko.inputs import signor

resources = signor()
```

After the first successful load, the cached release is available offline. Set
`NEKO_CACHE_DIR` to choose the cache root. Use `normalize_entities=False` only
when raw identifiers such as `SIGNOR-C1` are explicitly needed. Preloaded
DataFrames can also be supplied as `entity_dictionaries`.

ChEBI accessions are preserved as canonical non-protein identifiers rather
than being sent to UniProt. NeKo lazily caches the official compressed
`compounds.tsv.gz` table and extracts only the names required by the current
resource. Name enrichment is best-effort; offline or failed downloads fall
back to displaying the accession without blocking network construction.
[ChEBI data](https://www.ebi.ac.uk/chebi/) are provided by EMBL-EBI under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## PhosphoSitePlus identifiers

`neko.inputs.phosphosite()` preserves phosphorylation sites in the native
`GENE_RESIDUE` form, for example `MAP3K4_T1494`. Serine, threonine, and
tyrosine sites are recognized before protein identifier translation, so site
nodes are never submitted to UniProt. When a resource uses gene symbols as
edge identifiers, NeKo also keeps those symbols as the graph identifiers to
ensure kinase-to-site and site-to-protein paths remain searchable.
