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

---

## Class reference

::: neko._methods.enrichment_methods.Connections
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"
