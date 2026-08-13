# Strategy API

Connection strategies are exposed as methods on `Network`; their implementations
live in `neko.core.strategies`. Most users should call the bound methods so
history snapshots are recorded automatically.

See [Choosing a connection strategy](../strategies/index.md) for topology
diagrams, biological interpretation, and a decision table.

```python
net.connect_nodes(only_signed=True)
net.complete_connection(
    maxlen=2,
    path_policy="all_shortest",
    reuse_policy="discovered_paths",
    only_signed=True,
)
```

---

::: neko.core.strategies
    options:
      show_source: true
      show_root_heading: true
      heading_level: 2
      filters:
        - "!^_"
