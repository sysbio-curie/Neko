# History and network states

History is implemented directly by `Network`; there is no separate
`NetworkHistory` class. Mutating decorated methods save deep node/edge
snapshots as `NetworkState` objects. Undo, redo, checkout, and subsequent
mutations form a branching tree.

For structured topology explanations, provenance, pruning gaps, and JSON
export/recovery, see [Interaction curation](../interaction-curation.md).

```python
net.connect_nodes()
states = net.list_states()
net.undo()
net.redo()
html = net.history_html()
```

For routine inspection, avoid materializing complete history:

```python
status = net.history_status()
count = net.state_count
page = net.list_states(offset=0, limit=25, detail="summary")
details = net.get_state(page[0]["id"]) if page else None
```

`history_status()` returns `pruned_state_count` and
`working_state_matches_current`. It compares the working tables with the current
snapshot, including nested provenance, without traversing older snapshots.
It costs work proportional to the current network, not the retained history.
`state_count` is the number of retained states.

`list_states(*, offset=0, limit=None, detail="full")` slices IDs before copying
metadata. The default result is unchanged. Summary entries retain `id`,
`parent_ids`, `children_ids` and `current`; their metadata contains `method`,
`label`, `ancestry_gaps` and `transition`. The summary transition includes
`explanation_status` and `counts`, with `nodes`/`edges` counts for `added`,
`removed` and `changed`. Evidence, parameters and inclusion explanations are
available through `get_state`. Summary entries omit the repeated
`pruned_state_ids` list; use `history_status` for the count and `history_dict`
for exact IDs. All returned mutable values are independent copies.
For legacy snapshots without a recorded table delta, that table's counts are
`None`, rather than invented zeroes.

Offsets and limits must be non-negative integers (booleans are rejected);
`limit=None` selects the remaining states, and zero returns an empty list.
Unknown detail modes and invalid pagination values raise `ValueError`.
Pages use creation order; they are not frozen cursors across concurrent edits.

MCP summaries can use the page's parent/child IDs directly, together with
`state_count` and `history_status`. Calling `history_graph()` or `history_dict()`
as well would reintroduce full-history work.

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
        - state_count
        - root_state_id
        - history_graph
        - history_digraph
        - history_html
        - export_history
        - from_history
        - history_dict
        - history_status
        - get_state
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
