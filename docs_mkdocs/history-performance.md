# History performance and automatic network growth

This optimization keeps the existing `NetworkState` graph, schema version 1,
complete evidence, and independent snapshots. It does not prune history or
change connection policies.

## Changes

- `history_status()` compares only the working tables and current snapshot.
  `state_count` and paginated `list_states(detail="summary")` let MCP report
  history without materializing the whole history or a metadata-rich graph.
- Mutation detection reuses pre-operation fingerprints and avoids normalizing
  already normalized table records twice. Working-table fingerprints are not
  cached across operations, so direct nested evidence changes remain
  detectable. Immutable snapshots keep a SHA-256 digest of their tables, so the
  current snapshot is not serialized again by every later operation.
- Inclusion explanations use an exact-key index instead of scanning all earlier
  explanations. Actual distinct paths and operation ancestry remain recorded.
- Path cleanup first checks directed pairs. Unique pairs need no provenance
  serialization; duplicate pairs still receive exact full-record comparison.
- Single-row edge groups bypass parallel-evidence filtering while retaining
  effect/reference normalization and independent provenance. Conflicting manual
  effects and parallel resource evidence retain the existing merge behavior.
- Exporters capture independent working tables and one history payload rather
  than deep-copying the whole network, database resources and connection indexes.
  All Boolean variants reuse that capture. Compact JSON is written state by
  state; optional pretty printing remains available. Failed JSON serialization
  does not replace an existing destination.

## Paired measurements

Measured locally on 2026-10-20 against committed candidate
`ba87f7631dd6420b5269852daceb303f773e8963`, with the same Python 3.14 environment.
The optimized candidate is the worktree containing this report.
[Raw samples and output hashes](measurements/history-performance.json) are
included for review.

The resource is a synthetic directed ring with edges at offsets 1 and 5,
two edges per node. Connection starts with all resource nodes for
`connect_nodes`, and every tenth node for the other strategies. Completion uses
`all_shortest`, `discovered_paths`, and `maxlen=10`; radial growth uses three
outward levels; upstream growth uses depth/rank 3. These strategies select
different final networks, so their times are not a ranking of equivalent work.

Provider discovery and identifier mapping are offline. Construction of the seed
network is outside the timed strategy operation. Fresh processes use
`PYTHONHASHSEED=0` to make set-based frontier ordering reproducible. Medians use
five samples for the 100-node fixture and three for the 200-node fixture.
Wall-clock measurements on a shared host are indicative, not performance
guarantees or timing thresholds for CI.

| Resource nodes | Strategy, history enabled | Before (ms) | After (ms) | Reduction |
|---|---|---:|---:|---:|
| 100 | Connect existing nodes | 214.66 | 156.52 | 27.1% |
| 100 | Complete connection | 80.18 | 59.39 | 25.9% |
| 100 | Radial extension | 1,293.71 | 381.66 | 70.5% |
| 100 | Upstream extension | 97.99 | 66.43 | 32.2% |
| 200 | Connect existing nodes | 670.48 | 530.64 | 20.9% |
| 200 | Complete connection | 290.52 | 162.16 | 44.2% |
| 200 | Radial extension | 4,745.05 | 1,147.20 | 75.8% |
| 200 | Upstream extension | 190.68 | 126.67 | 33.6% |

For every paired strategy/size combination, both the complete edge-record hash
and the full history-payload hash match before and after optimization. This
checks effects, evidence, state relationships, effective parameters, inclusion
explanations and ordering, not just edge counts.

The raw benchmark also measures history-disabled execution as a diagnostic.
Do not subtract it from history-enabled time to estimate capture overhead:
the existing history-disabled implementation synchronizes object views at
nested calls, so the two modes have different execution paths.

## History inspection and exports

These fixtures import the full ring from SIF and apply label changes.
Measured label mutations themselves create additional states: the 100/200/500
node fixtures retain 26/24/56 states respectively.

| Nodes / edges | Operation | Before (ms) | After (ms) |
|---|---|---:|---:|
| 100 / 200 | Label mutation | 18.69 | 28.05 |
| 200 / 400 | Label mutation | 33.56 | 27.81 |
| 500 / 1,000 | Label mutation | 83.52 | 67.62 |
| 100 / 200 | SIF export with companion | 202.63 | 77.00 |
| 200 / 400 | SIF export with companion | 425.79 | 112.70 |
| 500 / 1,000 | SIF export with companion | 1,776.30 | 541.61 |

The small-fixture label timing regressed in this run; these measurements do not
establish a universal mutation speedup. Full `history_dict()` also remains
expensive: 393.30 ms before versus 396.53 ms after on the large fixture.
Its purpose is still to materialize every retained state.

On that same optimized large fixture, `history_status()` takes 11.75 ms and a
five-state summary page takes 0.018 ms. MCP needs to use these new APIs, including
page parent/child IDs, instead of also calling `history_dict()` or
`history_graph()` for each summary.

The large SIF companion shrinks from 35,442,561 to 15,811,104 bytes through compact
formatting. No snapshots or evidence fields are removed. Larger histories can
still exceed a downstream artifact-size limit.

## Validation and remaining limits

Regression tests exercise six representative strategies at two sizes. Each
outer operation creates one snapshot and makes at most six full-table
fingerprint calls (four once the current snapshot's digests are cached),
independent of how many nested additions occur. No full-history payload is
constructed inside these strategy calls. Tests also verify exact duplicate
handling, retained distinct paths, bounded page copying, nested evidence
mutation detection, atomic edits, exporter independence and history recovery.

Validation of this candidate:

- 257 offline library tests passed; 10 live/provider tests were deselected.
- The existing MCP real-NeKo Wave 2 workflow passed against this source checkout
  in Python 3.13, including Boolean/ODE handoff, provenance binding and recovery.
  Its import-time provider discovery required network access outside the sandbox.
- Wheel and sdist builds and strict Twine checks passed. A separate extracted
  wheel smoke test verified version 1.10.2, the summary APIs, automatic growth
  and SIF companion recovery.

These changes do not make growth cost independent of network size. Pair checks,
consolidation, path search and final snapshot capture still inspect network
data. Recording many selected paths or large evidence records also has a real
cost proportional to that information. The existing snapshot representation
still duplicates provenance across retained states; immutable sharing and a
versioned deduplicated storage format are separate future work.

Reproduce with the same interpreter in sequential processes:

```sh
python tests/benchmark_history_performance.py \
  --source /path/to/baseline --nodes 100 --steps 20 --output /tmp/before.json
python tests/benchmark_history_performance.py \
  --source /path/to/candidate --nodes 100 --steps 20 --output /tmp/after.json
```

Use `--nodes 200 --repeats 3` for the larger growth fixture and
`--nodes 500 --steps 50 --skip-strategies` for the large history/export fixture.
The script fixes its random hash seed before importing NeKo.
