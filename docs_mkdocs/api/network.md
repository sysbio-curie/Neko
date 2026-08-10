# Network

The `Network` class is the central object in NeKo. It holds a directed graph of biological nodes (genes, proteins, complexes) and edges (interactions), and exposes methods for expanding, connecting, querying, and exporting those graphs.

## Import

```python
from neko.core.network import Network
```

## Quick example

```python
from neko.core.network import Network
from neko.inputs import Universe

resources = Universe("omnipath")

net = Network(["EGFR", "KRAS", "MYC"], resources=resources.interactions)
net.connect_nodes()
print(net.nodes)
print(net.edges)
```

## Complete a seed network

`complete_connection` attempts both directed orientations for every original
seed pair. Path selection and reuse are explicit:

```python
net.complete_connection(
    maxlen=2,
    path_policy="all_shortest",
    reuse_policy="induced_subgraph",
    only_signed=True,
    consensus=False,
)
```

Path policies are `one_shortest`, `all_shortest`, and `all_bounded`. Reuse
policies are `none`, `discovered_paths`, and `induced_subgraph`. A finite
positive `maxlen` is mandatory. The legacy `algorithm`, `minimal`, and
`connect_with_bias` parameters remain temporarily available and emit a
migration warning with the equivalent explicit call.

If neither old nor new selectors are supplied, the transition release keeps
the former effective default: `all_bounded + discovered_paths`. An explicitly
disabled `minimal` flag maps to `none` unless bias is enabled; either biased
legacy combination maps to `induced_subgraph`.

See [Choosing a connection strategy](../strategies/index.md) for the full
policy matrix, fictitious topology diagrams, and the biological implications
of all public connection strategies.

## Connect to a GO term

```python
net.connect_genes_to_phenotype(
    id_accession="GO:0062043",
    only_signed=True,
    compress=True,
    maxlen=1,
)
```

The accession is sufficient: NeKo obtains the canonical term label from GO.
Exact-term human annotations are used by default. Use
`include_descendants=True` to include genes annotated to more specific GO
terms, or change `taxon_id` for another organism.

With `compress=True`, connected GO-associated genes are replaced by one node
named from the canonical GO term label. If collapsing those genes produces both
an activating and an inhibiting interaction between the same two nodes, NeKo
retains the conflicting evidence as one `bimodal` interaction. References and
interaction types from the contributing edges are preserved.

Parallel regulatory edges elsewhere in a network follow the same rule: an
`A stimulation B` edge together with an `A inhibition B` edge is represented
as `A bimodal B`. Complex formation remains separate because it does not encode
a regulatory sign.

---

## Class reference

::: neko.core.network.Network
    options:
      members:
        - __init__
        - add_node
        - add_edge
        - remove_node
        - connect_nodes
        - connect_subgroup
        - connect_component
        - connect_to_upstream_nodes
        - connect_genes_to_phenotype
        - complete_connection
        - is_connected
        - convert_edgelist_into_genesymbol
      show_source: true
      show_root_heading: true
      heading_level: 3
