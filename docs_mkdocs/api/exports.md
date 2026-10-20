# Exports

`Exports` converts a `Network` object into various file formats used by modelling tools such as MaBoSS, GINsim, and more.

## Import

```python
from neko._outputs.exports import Exports
```

## Supported formats

| Format | Method | Tool |
|---|---|---|
| BNet (Boolean Network) | `export_bnet()` | MaBoSS, PyBoolNet |
| SIF (Simple Interaction Format) | `export_sif()` | Cytoscape |

## Quick example

```python
from neko._outputs.exports import Exports

exporter = Exports(net)

# Export to BNet format for MaBoSS
exporter.export_bnet("my_model.bnet")

# Export to SIF for Cytoscape
exporter.export_sif("my_network.sif")
```

`export_bnet("my_model.bnet")` writes numbered files such as
`my_model_1.bnet`, because each bimodal edge can produce stimulation and
inhibition variants. Opposite parallel regulatory edges are normalized to one
bimodal edge before export. Pass `n=` to cap the number of variants generated.

SIGNOR group/context labels and other identifiers containing spaces, `/`, `-`,
`#`, or `:` are sanitized for BNet output. Export raises `ValueError` if two
different labels would collapse to the same sanitized identifier. Custom nodes,
including compressed GO phenotypes, are resolved using the network's node table
rather than being treated as protein identifiers. Null endpoints and endpoints
missing from that node table are rejected with a descriptive error.

Each output also receives a `<filename>.neko.json` companion preserving source
state identity, canonical-to-exported identifier mappings, provenance, and retained
history. SIF additionally embeds provenance comments. Keep the JSON companion for
formats such as BNET that cannot carry this metadata. See
[Interaction curation and persistence](../interaction-curation.md).

`Exports(net)` captures independent working tables and one history payload when
constructed, without copying resource tables or connection indexes. Later edits,
checkout or pruning of `net` do not change that export capture. All Boolean
variants reuse the captured history.

Companions and `net.export_history(path)` use compact JSON by default, retaining
the same schema and every snapshot/evidence field. Use `Exports(net, pretty=True)`
or `net.export_history(path, pretty=True)` for indented JSON. Compact serialization
writes individual states without allocating one full-history JSON string.
Snapshot retention and downstream artifact-size limits still apply.

---

## Class reference

::: neko._outputs.exports.Exports
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"
