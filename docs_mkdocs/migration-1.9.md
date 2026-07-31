# Upgrading to NeKo 1.9

NeKo 1.9 is the final major feature release planned for the current backend.
It modernizes identifier translation, SIGNOR ingestion, Gene Ontology
queries, tissue mapping, and network export. The existing tutorial notebooks
and principal `Network` workflows remain valid.

Version 2.0.0 is reserved for the separately developed backend rewrite.
NeKo 1.9 intentionally does not restore compatibility wrappers for the
localized interfaces removed below.

## Install or upgrade

The distribution is named `nekomata`, while the import package remains
`neko`:

```bash
python -m pip install --upgrade nekomata==1.9.0
```

Verify the installed version:

```python
import importlib.metadata
import neko

assert importlib.metadata.version("nekomata") == "1.9.0"
assert neko.__version__ == "1.9.0"
```

## Identifier translation and caching

NeKo no longer uses PyPath for gene-symbol and UniProt translation. It first
uses a cached table of reviewed human UniProt entries and then falls back to
the official UniProt REST mapping service for identifiers not found locally.

The first translation may therefore require network access. Later lookups use
the local cache. Set `NEKO_CACHE_DIR` before starting Python to choose its
location:

```bash
export NEKO_CACHE_DIR=/path/to/neko-cache
```

The translation API is available directly when needed:

```python
from neko.inputs import identifier_mapping

uniprot = identifier_mapping.to_uniprot("TP53")
symbol = identifier_mapping.to_genesymbol("P04637")
```

Lightweight `pypath-common` utilities remain a dependency, but
`pypath-omnipath`, PyPath mapping, and the associated Paramiko dependency
chain have been removed.

## Gene Ontology interfaces

The following documented interfaces were removed:

- `neko._annotations.gene_ontology.fetch_nodes_from_url`
- `Ontology.modify_url_ontology`

Use the official GO API methods on `Ontology` instead:

```python
from neko._annotations.gene_ontology import Ontology

ontology = Ontology(taxon_id=9606)
term = ontology.get_term("GO:0062043")
genes = ontology.fetch_go_genes("GO:0062043")
markers = ontology.get_markers(id_accession="GO:0062043")
```

GO accessions are authoritative. Prefer `id_accession=` when connecting a
network to a phenotype. Free-text phenotype names resolve only when they are
registered as local aliases.

Exact-term associations are returned by default. To include associations
propagated from descendant terms:

```python
markers = ontology.get_markers(
    id_accession="GO:0062043",
    include_descendants=True,
)
```

To exclude automatic assertions:

```python
markers = ontology.get_markers(
    id_accession="GO:0062043",
    exclude_automatic_assertions=True,
)
```

Unknown accessions raise `GeneOntologyNotFoundError`. Transport, decoding,
pagination, and response-schema failures raise `GeneOntologyError`.

## SIGNOR normalization

`signor()` now downloads and validates the SIGNOR interaction table and its
entity dictionaries through NeKo's managed cache. Complexes, protein
families, phenotypes, and stimuli are normalized to readable typed
identifiers by default:

```python
from neko.inputs import signor

resources = signor()
```

Examples of normalized prefixes include:

- `COMPLEX:`
- `PROTEIN_FAMILY:`
- `PHENOTYPE:`
- `STIMULUS:`

If an existing workflow requires raw proprietary SIGNOR identifiers, disable
normalization explicitly:

```python
resources = signor(normalize_entities=False)
```

SIGNOR's directness field now describes evidence directness separately from
the direction of the graph edge. Code reading raw SIGNOR-derived columns
should account for the new `is_direct` field.

## Tissue mapping

Tissue mapping validates Human Protein Atlas annotations before classifying
genes. Cancer tissues use a managed cache of the HPA cancer dataset; other
tissues use HPA annotations obtained through OmniPath.

Service and schema failures raise `AnnotationServiceError` rather than
silently classifying genes as absent.

## BNet and SIF exports

Opposing or duplicate regulatory evidence is consolidated before export while
preserving references. Parent directories are created automatically for BNet
and SIF destinations.

BNet identifiers are sanitized for BoolNet compatibility. If two distinct
labels collapse to the same identifier, export now raises `ValueError`
instead of writing an ambiguous model. Rename the conflicting nodes before
exporting.

Each bimodal edge can generate stimulation and inhibition variants. Use `n`
to limit the number of files:

```python
from neko._outputs.exports import Exports

Exports(network).export_bnet("models/network.bnet", n=16)
```

## Compatibility summary

Most users following the bundled notebooks do not need code changes. Review
an existing workflow if it:

- imported either removed ontology interface;
- depended on raw SIGNOR identifiers;
- interpreted SIGNOR directness as graph direction;
- expected BNet export to tolerate sanitized-name collisions; or
- assumed that identifier translation always contacted PyPath.

See the [changelog](changelog.md) for the complete release summary.
