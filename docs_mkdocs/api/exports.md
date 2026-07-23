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

---

## Class reference

::: neko._outputs.exports.Exports
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"
