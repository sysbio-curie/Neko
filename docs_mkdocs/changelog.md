# Changelog

All notable changes to NeKo are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and version numbers generally follow [Semantic Versioning](https://semver.org/).
NeKo 1.9.0 is a transitional release with a small number of explicitly
documented breaking changes; version 2.0.0 is reserved for the separately
developed backend rewrite.

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
- Branching **NetworkHistory** with automatic state snapshots and HTML/SVG rendering
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
