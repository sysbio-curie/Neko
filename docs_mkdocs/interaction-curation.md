# Interaction curation and reproducible history

NeKo owns interaction validation, provenance, and state transitions. An MCP
wrapper should pass authorized requests to `Network.edit_interactions` and
return its actual result. It should not edit `net.edges` directly or reproduce
NeKo's sign, identifier, consolidation, or history logic.

## Public editing contract

```python
from neko.core.network import InteractionEditError

result = net.edit_interactions(
    edits=[{
        "operation": "add",
        "source": "P31749",
        "target": "Q15796",
        "effect": "inhibition",
        "rationale": "Explain why the cited interaction belongs in this model.",
        "evidence": {
            "summary": "Describe the specific observation supporting this effect.",
            "pmids": [],
            "citations": ["Supply the actual supporting publication here"],
        },
    }],
    expected_state_id=net.current_state_id,
)
```

The example's citation is a placeholder, not scientific evidence. Substitute
the actual supporting publication when using the API.

Each edit has exactly the six keys shown. `operation` is `add` or `set_effect`.
The public `InteractionEdit`, `InteractionEvidence`, `InteractionChange`, and
`InteractionEditResult` TypedDict definitions are in `neko.core.curation`.
An `add` requires an absent directed pair; `set_effect` requires exactly one
existing interaction for that pair. A batch may address each pair once. All
requests are validated before any changes are published, including the new
snapshot. A rejected batch changes neither topology nor history.

Endpoints resolve against the current node table's canonical identifiers and
labels. An identifier matching multiple canonical nodes is ambiguous. No
provider lookup, implicit node creation, or automatic reconnection occurs.
Manual interactions need not occur in the selected database. Custom-node
admission remains a separate operation.

Allowed effects are `stimulation`, `inhibition`, `bimodal`, `form complex`, and
`undefined`. Rationale and evidence summary must be nonempty strings. Evidence
requires at least one PMID or other citation. PMIDs are positive decimal strings
without a `PMID:` prefix; absent PMIDs are `[]`. The library checks structure,
not publication existence, scientific validity, or researcher authorization.

The JSON-compatible result is:

```json
{
  "changed": true,
  "previous_state_id": 12,
  "state_id": 13,
  "changes": [
    {
      "operation": "add",
      "source": "P31749",
      "target": "Q15796",
      "before": null,
      "after": {"source": "P31749", "target": "Q15796", "Effect": "inhibition"}
    }
  ]
}
```

The `after` record above is abbreviated: actual records also include `Type`,
`References`, and the complete `Provenance`. `before` contains the full previous
record for an effect change. One changed batch creates one state. Setting an
already-identical effect returns `changed: false`, the unchanged state ID, and
an empty `changes` list; it does not append evidence. Mixed batches report only
actual changes.

`expected_state_id` is an optional optimistic guard against an unexpected
current state. State IDs are local to a network history, not global identifiers.
The MCP layer must also bind requests to its network/session identity and
authorization. Any changed result invalidates handoff approvals associated with
the previous state; NeKo does not grant approvals or export automatically.

### Errors

`InteractionEditError` is a `ValueError` with `code`, `edit_index` (zero-based, or
`None` for a batch-level error), and JSON-compatible `details` attributes.

| Code | Meaning |
| --- | --- |
| `invalid_request` | Missing/extra keys, invalid operation, empty rationale, repeated pair, or invalid request shape |
| `invalid_effect` | Effect is outside the supported vocabulary |
| `invalid_evidence` | Evidence is missing or malformed |
| `endpoint_not_found` | Endpoint is not in the current network |
| `ambiguous_endpoint` | Identifier resolves to multiple existing nodes |
| `interaction_exists` | Addition would overwrite an existing directed interaction |
| `interaction_not_found` | Effect-change target is absent |
| `ambiguous_interaction` | Multiple interactions share the directed pair |
| `stale_state` | Current state differs from `expected_state_id` |
| `history_disabled` | Automatic history is suspended or disabled |

No mutation success should be reported after an exception. The older low-level
`add_edge` remains available for construction compatibility, but it is not the
manual-curation interface. A rejected implicit endpoint now raises `ValueError`
and rolls back that low-level call's endpoint additions.

## Edge provenance

Existing edge columns keep their meanings. The new `Provenance` column contains
a dictionary, also available as `Edge.metadata["Provenance"]`:

| Field | Meaning |
| --- | --- |
| `origin` | Original introduction through `resource`, `manual`, `import`, or `unknown` |
| `original` | Original `origin`, `database_name`, `effect`, and `references` |
| `originals` | Present when consolidation retains multiple original evidence records |
| `edits` | Ordered manual operations with `origin`, `state_id`, `before_effect`, `after_effect`, `rationale`, and `evidence` |

Changing a database interaction's effect keeps its resource origin and original
references; the change itself is explicitly manual. Adding an absent interaction
always records manual origin, even if the same pair happens to exist in a database.
Low-level supplied rows of unknown provenance are not relabelled as database
evidence. Unannotated SIF input has import origin.

Available manual PMIDs are retained in the edit evidence; for a new manual edge
they also populate `References`. An effect change preserves `References` rather
than replacing it with edit citations. Parallel database evidence is retained
in `originals` without changing the established working-sign selection policy.
Subsequent consolidation retains curated effects and their explanations. Two
conflicting curated effects cannot be silently consolidated.

Because provenance is structured, consumers should select the legacy columns
explicitly when using hash-based DataFrame operations or five-column tuple
comparisons. `consolidate_edges` handles provenance itself. Mutating returned
result dictionaries, copied networks, or restored working frames does not mutate
historical snapshots.

## Database identity

`Network(..., database_name="SIGNOR")` names a supplied resource table. Built-in
loaders carry their names into `Universe`, its interaction DataFrame attributes,
and the constructed network. `net.database_name` is read-only and captured in
every state. An explicit name conflicting with a known loader name is rejected.
Unnamed input remains `None`; a selected database does not establish every edge's
origin. A later change to a loader/session default does not relabel the network.

Only the database name is recorded here. Database version discovery, acquisition
times, dataset hashes, and database snapshots are outside this feature.

## Existing history and explanations

The existing `save_state`, `checkout`, `undo`, `redo`, `compare_states`,
`list_states`, `history_graph`, `history_digraph`, and `history_html` APIs remain.
Checkout restores node/edge objects, structured provenance, and initial seeds as
well as DataFrames. Comparison retains its original summary keys and adds exact
`node_changes` and `edge_changes` records with `added`, `removed`, and `changed`
lists; changed records have `before` and `after` values.

`get_state(state_id=None)` reads one exact state ID (default: current) without
checkout. It returns `id`, `current`, `parent_ids`, `children_ids`, `metadata`, and
`nodes`/`edges` tables, each represented by `columns` and `records` lists.
`history_dict()` returns the same independent JSON-compatible dictionary written
by `export_history`; callers need not access private fields or parse display text.
For lightweight inspection, use `history_status()` for pruning and working-state
flags, `state_count` for the retained count, and
`list_states(offset=0, limit=25, detail="summary")` for a page without full evidence.
See [History API](api/history.md) for the exact summary contract.

State metadata adds `database_name`, `initial_nodes`, and `transition`:

- `parent_state_id`: the actual parent when the transition was created.
- `nodes` and `edges`: exact snapshot differences, including evidence changes.
- `operations`: construction calls with effective parameters, including defaults
  and resolved path/reuse policies. Nested calls retain their actual parameters.
- `inclusion_reasons`: selected paths, shortest-path unions identified by seed
  pair, induced-subgraph membership, upstream cascade edges, ontology association
  or compression, and manual edit evidence. `operation_indices` references entries
  in `operations`, retaining nested context without repeating full parameters.
- `explanation_status` and `limitations`: distinguish recorded, partial, and
  unavailable explanations.

The all-shortest algorithm selects an edge union. Its records preserve that union
and seed pair; they do not invent a unique path. All-bounded selections can record
several actual supporting paths for the same edge. Existing upstream-cascade
construction does not apply its advertised sign filter; its effective record
explicitly reports `sign_filter_applied: false` rather than changing that existing
scientific behavior in this feature.

Display arguments remain concise. Authoritative `operations`, snapshots, and
manual explanations are not truncated. The HTML state graph links to expandable
structured details. Metadata strings are escaped for HTML rendering.

Retention behavior is unchanged. `ancestry_gaps` and `pruned_state_ids` report
removed states, while `transition.parent_state_id` retains the original parent
even if the visible graph reconnects retained ancestors. Imported or untracked
changes report missing explanations explicitly. Direct DataFrame writes cannot
be assigned scientific explanations after the fact.

There is intentionally no single flag promising complete reproducibility. A
summary must inspect `pruned_state_ids`, state-level `ancestry_gaps`, each
`transition.explanation_status`/`limitations`, and `working_state_matches_current`.
Older serialized NeKo objects retain their existing snapshots and explicitly
report unavailable explanations and unknown prior retention gaps on loading.

## Persistence and graph exports

```python
net.export_history("network.history.json")
restored = Network.from_history("network.history.json")
restored.checkout(saved_state_id)

# Explicitly supply resources to perform further database-driven construction.
restored = Network.from_history("network.history.json", resources=resource_table)
```

JSON is compact by default; `export_history(path, pretty=True)` and
`Exports(net, pretty=True)` enable indentation without changing the schema.

The versioned `neko-history` JSON format (`schema_version: 1`) stores retained
snapshots, metadata, parent/child relationships, state counter, current/root IDs,
retention settings, pruned-state IDs, database name, and the working network.
`working_state_matches_current` reports whether unrecorded working changes exist.
They are preserved separately rather than falsely attributed to the current state.
Loading validates graph relationships and endpoint consistency. It does not
download database data. Offline restoration permits navigation, export, and
explicit manual curation; connection construction requires supplied resources.

This restores stored topology and evidence. It does not prove that construction
will replay identically against a different database with the same name.

Every SIF or BNET output has a `<filename>.neko.json` companion (`neko-export`,
schema version 1) with the source history, state ID, identifier mappings, output
SHA-256, and export interaction records. BNET records identify the chosen bimodal
variant; undefined effects are reported as omitted. The hash identifies the graph
artifact, not a database dataset. Load the companion using `Network.from_history`
to recover the original NeKo network, not the Boolean projection.

SIF also embeds NeKo provenance comments which its importer preserves. A plain
SIF import creates a new history; preserved edit IDs refer to the source history.
Use the JSON companion to recover original branch relationships and exact
canonical identifiers. Other SIF consumers may ignore comments. BNET cannot
encode this provenance itself, so retain its companion when transferring it.

## MCP integration boundary

Expose a focused `edit_interactions` tool accepting this schema, pass the expected
state, and return NeKo's resulting canonical identifiers and actual changes.
Translate `InteractionEditError` without retrying as an implicit node addition.
Keep researcher authorization, session revision/approval checks, artifact
registration, and result pagination in the MCP layer. Use existing history tools
to inspect `transition` and expose a reference to `export_history` for full details.
No parallel history implementation, database-sign reimplementation, or additional
publication authorization is needed in NeKo.
