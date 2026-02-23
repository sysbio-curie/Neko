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

resources = Universe()
resources.build()

net = Network(["EGFR", "KRAS", "MYC"], resources=resources.interactions)
net.connect_nodes()
print(net.nodes)
print(net.edges)
```

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
