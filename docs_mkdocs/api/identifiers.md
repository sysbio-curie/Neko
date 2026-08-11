# Identifier mapping

NeKo lazily loads a cached reviewed-human UniProt table for gene-symbol and
accession translation. Unrecognized identifiers can use a bounded live UniProt
fallback; successful and failed lookups are memoized for the process.

```python
from neko.inputs.identifier_mapping import to_genesymbol, to_uniprot

accession = to_uniprot("EGFR")
symbol = to_genesymbol("P00533")
```

Set `NEKO_CACHE_DIR` to choose the cache root. `refresh_cache()` explicitly
refreshes the reviewed-human table.

::: neko.inputs.identifier_mapping.looks_like_uniprot_accession

::: neko.inputs.identifier_mapping.to_uniprot

::: neko.inputs.identifier_mapping.to_genesymbol

::: neko.inputs.identifier_mapping.refresh_cache

## Batched translator

`IDTranslator` is the separate batched translation interface used by legacy
workflows.

::: neko.inputs.db_translator.IDTranslator
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"
