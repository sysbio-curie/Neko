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
