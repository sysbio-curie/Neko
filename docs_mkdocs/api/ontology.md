# Ontology

The `Ontology` class provides Gene Ontology (GO) utilities: fetching gene sets from GO accessions and mapping phenotypic categories to network nodes.

## Import

```python
from neko._annotations.gene_ontology import Ontology
```

## Quick example

```python
from neko._annotations.gene_ontology import Ontology

onto = Ontology()

# Fetch all genes annotated to a specific GO term
genes = onto.fetch_genes_from_go_id("GO:0007165")  # signal transduction

# Map nodes in a network to tissue expression data
onto.add_tissue_expression(net, tissue="liver")
```

---

## Helper function

::: neko._annotations.gene_ontology.fetch_nodes_from_url
    options:
      show_root_heading: true
      heading_level: 3

---

## Class reference

::: neko._annotations.gene_ontology.Ontology
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"
