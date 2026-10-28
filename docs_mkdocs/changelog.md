# Changelog

All notable changes to NeKo are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and version numbers generally follow [Semantic Versioning](https://semver.org/).
NeKo 1.9.0 is a transitional release with a small number of explicitly
documented breaking changes; version 2.0.0 is reserved for the separately
developed backend rewrite.

---

## [1.10.3] – 2026-10-28

### Added

- `Network.seed_mapping` records what happened to each requested seed: its
  database symbol (aliases such as `FAK` map to `PTK2`) or `None` when it was
  dropped. The initial history state stores the same mapping.
- `NetworkVisualizer.render()` accepts `format=` (`"pdf"` by default; also
  `"svg"`, `"png"` or any other Graphviz format).

### Fixed

- `Network(genes)` without `resources` loads OmniPath, as documented, instead of
  failing with `TypeError: 'NoneType' object is not subscriptable`.
- Custom interaction tables: text labels in an `effect` column are parsed by
  meaning. Before this fix, any non-empty text, including `"inhibition"`, was read as a
  stimulation, and missing effects were read as inhibitions. Unrecognised effect
  values now raise `ValueError` that names them, and missing values are
  `undefined`. `form complex` sets `form_complex`.
- Boolean flag columns (`is_stimulation`, `is_inhibition`, `is_directed`,
  `form_complex`) accept `TRUE`/`FALSE`, `t`/`f` and `y`/`n` in any case, and
  raise `ValueError` for other text instead of reading it as `True`.
- Seeds and nodes missing from the database are reported with a `UserWarning`
  instead of a printed message, so they can be caught, logged or turned into
  errors.
- Documentation: corrected the `render()` example and the default-resource note
  in the network-building notebook, and documented the accepted values for
  custom interaction tables.

---

## [1.10.2] – 2026-10-21

### Added

- Atomic, validated interaction additions and effect changes between existing
  endpoints, with manual origin, evidence, PMIDs and edit explanations.
- Structured topology explanations in the existing network history, JSON
  recovery, and provenance companions for SIF and Boolean exports.
- Selected database identity on networks and regression coverage for mutation
  validation, provenance preservation and history consistency.
- Lightweight history status, retained-state count and paginated summary access.

### Changed

- Avoid repeated full-table serialization during path deduplication and repeated
  inclusion-reason scans during automatic network growth.
- Reuse mutation fingerprints and normalize unique edge groups without parallel
  evidence merging overhead.
- Capture export history once without copying database resources; write compact
  JSON companions by default with optional pretty printing.

### Fixed

- Array-like arguments (NumPy arrays, pandas Series/Index, frozensets) no longer
  make recorded operations fail; unsupported values are recorded by type and
  repr instead of rejecting the call.
- A failure while recording history no longer leaves a stale transition context
  for a later `save_state()`.
- SIF import treats exported `None`/`nan` reference comments as missing and uses
  the `SIF file` fallback instead of storing them as references.
- Deep-copied history tables keep nested provenance when index labels repeat.
- Networks pickled before history tracking can grow again; their next change
  becomes the root state.
- `describe_history()` and `history_digraph()` hide internal explanation
  metadata and show user-supplied metadata again.
- History and export companion JSON files follow the umask instead of 0600.
- Multi-variant BNet exports write the history once to
  `<stem>.bnet.history.json`; their companions reference it (`neko-export`
  schema version 2) instead of each embedding a full copy. Single-model exports
  keep the self-contained schema version 1 companion.

---

## [1.10.1] – 2026-10-01

### Fixed

- SIF import keeps each edge's `# Reference PMID:` comment as its reference
  across export/import round trips. Unannotated interactions keep the
  `SIF file` fallback, and a comment is never carried over to another edge.

---

## [1.10.0] – 2026-08-13

NeKo 1.10.0 improves the performance, reproducibility, and configurability of
network completion while preserving compatibility with legacy connection
parameters through explicit migration warnings.

### Added

- Explicit `path_policy` values (`one_shortest`, `all_shortest`, and
  `all_bounded`) for `complete_connection`.
- Explicit `reuse_policy` values (`none`, `discovered_paths`, and
  `induced_subgraph`).
- Shortest-path predecessor-DAG selection for the union of all shortest paths.
- Visible migration warnings with exact replacements for legacy connection
  parameters.
- Expanded API documentation and regression coverage for connection policies.

### Changed

- `complete_connection` now checks the two directed orientations sequentially,
  allowing newly exposed paths to prevent redundant resource searches.
- New completion policies require a finite positive `maxlen`.
- Resource edge lookup, path insertion, cascade insertion, and induced closure
  use indexed and batched mutation.
- Neighbor ordering is stable for reproducible unweighted shortest-path choice.
- Ontology phenotype and pathway-recreation tutorials are reproducible with the
  current APIs.

### Fixed

- Restored the semantic distinction between independent searches and reuse of
  discovered paths.
- Working `Effect="undefined"` edges are no longer treated as signed paths.
- Removed deduplication of stale DataFrame objects after path mutation.
- History metadata no longer attempts identifier translation for policy strings.
- Node renaming changes display labels without rewriting canonical identifiers.

### Deprecated

- `algorithm`, `minimal`, and `connect_with_bias` in `complete_connection`; use
  `path_policy` and `reuse_policy`.

### Validation and compatibility

- The refactor passed 148 deterministic tests, including every legacy
  BFS/DFS × `minimal` × `connect_with_bias` mapping and the corresponding
  explicit-policy topology.
- On the pinned 18-gene SIGNOR benchmark, performance-only changes preserved
  100% of nodes, signed directed edges, initial seeds, and seed-pair directed
  reachability for `connect_nodes`, `all_bounded + discovered_paths`,
  `connect_subgroup`, `connect_component(mode="OUT")`, upstream connection,
  and Atopo-complete construction.
- `all_bounded + discovered_paths` improved from 113.9 seconds to 5.3 seconds
  (21.6×), while retaining its 111-node/771-edge topology exactly.
- Historical BFS selected equal-length routes through unordered set iteration.
  The deterministic `one_shortest` policy can therefore choose a different
  but equally short bridge topology; it retained 100% of initial seeds and
  directed seed-pair reachability in the reference benchmark.
- For users who need robustness across all equal-length alternatives,
  `all_shortest + none` contained 100% of the historical BFS topology, while
  `all_shortest + discovered_paths` retained 97.3% of historical nodes and
  95.7% of historical signed edges. Both retained 100% of seeds and directed
  reachability and remained 6.6–8.0× faster than historical BFS.

---

## [1.9.0] – 2026-07-31

NeKo 1.9.0 is the final major feature release planned for the current
backend. Existing tutorial notebooks and the principal `Network` workflows
remain supported. See [Upgrading to NeKo 1.9](migration-1.9.md) for migration
details.

### Added

- Cache-backed reviewed-human UniProt identifier mapping with bounded,
  retried downloads and a UniProt REST fallback.
- A shared cache location controlled through `NEKO_CACHE_DIR`.
- Validated SIGNOR database caching and entity dictionaries for complexes,
  protein families, phenotypes, and stimuli.
- SIGNOR ChEBI display-name resolution backed by the official ChEBI compound
  table.
- Official Gene Ontology API integration with validated pagination, taxon
  filtering, descendant controls, and optional removal of automatic
  assertions.
- Structured `GOTerm` and `GOGene` records.
- Deterministic regression coverage for identifier mapping, SIGNOR, Gene
  Ontology, tissue mapping, ChEBI, phosphosites, and export behavior.

### Changed

- Replaced PyPath-based identifier translation and removed the
  `pypath-omnipath`/Paramiko dependency chain; lightweight `pypath-common`
  utilities remain.
- SIGNOR entities are normalized by default to readable typed identifiers
  such as `COMPLEX:`, `PROTEIN_FAMILY:`, `PHENOTYPE:`, and `STIMULUS:`.
- SIGNOR directness evidence is stored separately from graph direction.
- Gene Ontology phenotype lookup now uses authoritative GO accessions and the
  official GO API instead of the former URL-based workflow.
- Tissue expression queries use validated Human Protein Atlas annotations,
  with managed caching for the HPA cancer dataset.
- Duplicate and opposing regulatory evidence is consolidated before export.
- BNet variant generation is lazy and can be bounded with `n=`.
- Runtime dependency constraints support pandas 2 and 3 and current Python
  releases.
- Package maturity metadata is now Beta.

### Fixed

- Preserved custom phenotype labels during identifier translation and export.
- Prevented ambiguous BNet models when distinct labels collide after
  identifier sanitization.
- Rejected null, empty, or unknown BNet endpoints with descriptive errors.
- Created parent directories automatically for BNet and SIF exports.
- Preserved references while consolidating bimodal regulatory evidence.
- Normalized phosphosite identifiers consistently.
- Retried and rejected invalid or truncated SIGNOR downloads.
- Corrected installed version discovery to query the `nekomata`
  distribution.

### Breaking changes

- Removed the documented
  `neko._annotations.gene_ontology.fetch_nodes_from_url` helper.
- Removed `Ontology.modify_url_ontology`; use
  `Ontology.get_term`, `Ontology.fetch_go_genes`, or
  `Ontology.get_markers` with a GO accession.
- `signor()` now normalizes proprietary SIGNOR entities by default. Pass
  `normalize_entities=False` when raw SIGNOR identifiers are required.
- BNet export now raises `ValueError` for identifiers that become ambiguous
  after sanitization instead of writing an invalid model.

No compatibility wrappers are provided for the removed ontology interfaces.

---

## [1.1.0] – 2025

### Added
- Branching network history with automatic `NetworkState` snapshots and HTML/SVG rendering
- `NetworkState` class for point-in-time network snapshots
- BFS / DFS graph traversal algorithms in `neko.core.algorithms`
- `connect_to_upstream_nodes` method
- Performance benchmarking scripts

### Changed
- Interaction lookup tables pre-processed for O(1) neighbour queries
- `connect_nodes` signature updated for clarity
- `pandas` pinned to `2.2.2` for stability

### Fixed
- Edge colouring bug in `NetworkVisualizer.vis_comparison`
- Handling of complex node names containing colons

---

## [1.0.0] – 2024

### Added
- Initial public release on PyPI as **nekomata**
- Core `Network` class with `add_node`, `add_edge`, `remove_node`
- `connect_nodes`, `connect_subgroup`, `connect_component`, `complete_connection`
- `connect_genes_to_phenotype` via Gene Ontology
- `NetworkVisualizer` with Graphviz and yFiles backends
- `Exports` class: BNet and SIF formats
- `Ontology` class for GO-term and tissue mapping
- 11 tutorial notebooks
- Sphinx documentation hosted on GitHub Pages

---

## [0.x] – Pre-release

Internal development versions at Institut Curie / Sysbio-Curie.
