# History and network states

History is implemented directly by `Network`; there is no separate
`NetworkHistory` class. Mutating decorated methods save deep node/edge
snapshots as `NetworkState` objects. Undo, redo, checkout, and subsequent
mutations form a branching tree.

```python
net.connect_nodes()
states = net.list_states()
net.undo()
net.redo()
html = net.history_html()
```

## Network history methods

::: neko.core.network.Network
    options:
      members:
        - save_state
        - set_max_history
        - set_history_tracking
        - suspend_history
        - list_states
        - checkout
        - restore_state
        - undo
        - redo
        - compare_states
        - describe_history
        - describe_states
        - current_state_id
        - root_state_id
        - history_graph
        - history_digraph
        - history_html
      show_source: true
      show_root_heading: true
      heading_level: 3

## Snapshot value object

::: neko.core.network_state.NetworkState
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"

## Rendering helpers

::: neko._visual.history.build_history_graph

::: neko._visual.history.history_digraph

::: neko._visual.history.history_html
