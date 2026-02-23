# Exports

`Exports` converts a `Network` object into various file formats used by modelling tools such as MaBoSS, GINsim, and CoBrexa.

## Import

```python
from neko._outputs.exports import Exports
```

## Supported formats

| Format | Method | Tool |
|---|---|---|
| BNet (Boolean Network) | `export_bnet()` | MaBoSS, PyBoolNet |
| SIF (Simple Interaction Format) | `export_sif()` | Cytoscape |
| GML | via `networkx` on `net.graph` | Various |
| GraphML | via `networkx` on `net.graph` | Gephi, yEd |

## Quick example

```python
from neko._outputs.exports import Exports

exporter = Exports(net)

# Export to BNet format for MaBoSS
exporter.export_bnet("my_model.bnet")

# Export to SIF for Cytoscape
exporter.export_sif("my_network.sif")
```

!!! tip "Exporting directly from `Network`"
    The convenience wrappers on `Network` call `Exports` internally, so in most cases you do not need to instantiate `Exports` yourself:

    ```python
    net.export_bnet("model.bnet")
    ```

---

## Class reference

::: neko._outputs.exports.Exports
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"
