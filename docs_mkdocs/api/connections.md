# Connections

`Connections` provides the search-and-connect algorithms that underlie `Network` expansion methods. It is initialised with an interaction database DataFrame and pre-processes lookup tables for fast neighbour queries.

You rarely need to instantiate `Connections` directly — it is used internally by `Network`. The documentation here is aimed at developers who want to extend NeKo with custom connection strategies.

## Import

```python
from neko._methods.enrichment_methods import Connections
```

## Quick example

```python
import pandas as pd
from neko._methods.enrichment_methods import Connections

db = pd.read_csv("my_interactions.csv")   # source, target, effect, ...
conn = Connections(db)

# Check if a direct path exists between two proteins
paths = conn.find_paths("EGFR", "AKT1", maxlen=3)
```

## Search semantics

The public `Network.complete_connection` API describes results through
`path_policy` rather than exposing traversal details:

- `one_shortest` uses BFS and selects one stable minimum-edge path.
- `all_shortest` uses a BFS predecessor DAG and selects the edge union of all
  minimum-edge paths.
- `all_bounded` uses bounded DFS and selects all simple paths through the
  cutoff.

All public completion policies require a finite positive cutoff. Low-level
`Connections.bfs` retains its legacy `force` behavior for specialized internal
use, but unbounded traversal is not an implicit network-construction policy.

`Connections` also stores indexed resource rows and signed adjacency maps.
Connection strategies should use those indexes and bulk network mutation rather
than scanning a DataFrame or calling `Network.add_edge` for every interaction.

---

## Class reference

::: neko._methods.enrichment_methods.Connections
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"
