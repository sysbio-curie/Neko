# Changelog

All notable changes to NeKo are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and the project adheres to [Semantic Versioning](https://semver.org/).

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
