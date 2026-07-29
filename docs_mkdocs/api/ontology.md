# Ontology

The `Ontology` class retrieves GO term metadata and gene associations from
the official Gene Ontology API. GO accessions are authoritative; phenotype
text is used only as a display label or a locally registered alias.

## Import

```python
from neko._annotations.gene_ontology import Ontology
```

## Quick example

```python
from neko._annotations.gene_ontology import Ontology

onto = Ontology(taxon_id=9606)

# Structured records retain both the gene symbol and source identifier.
genes = onto.fetch_go_genes("GO:0062043")
print([(gene.symbol, gene.gene_id) for gene in genes])

# The backward-compatible helper returns symbols only.
markers = onto.get_markers(id_accession="GO:0062043")
```

By default, only associations whose object is the requested GO term are
returned. Set `include_descendants=True` to include annotations propagated
from more specific terms. The default taxon is human (`NCBITaxon:9606`).
Numeric taxonomy IDs and complete `NCBITaxon:` CURIEs are both accepted.

Automatic assertions are included by default. Pass
`exclude_automatic_assertions=True` to remove `ECO:0000501` associations.
This and the exact-term filter are enforced locally because deployments of
the upstream API have not always applied their corresponding query flags.

HTTP, decoding, and response-schema failures raise `GeneOntologyError`.
Unknown accessions raise `GeneOntologyNotFoundError`; a valid term with no
matching genes returns an empty list.

---

## Class reference

::: neko._annotations.gene_ontology.Ontology
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"
