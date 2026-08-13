# Core objects and utilities

The working `Network.nodes` and `Network.edges` DataFrames are authoritative.
`Node` and `Edge` are value-object views that are synchronized explicitly by
network operations.

## Node

::: neko.core.node.Node
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"

## Edge

::: neko.core.edge.Edge
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"

## Public graph utilities

::: neko.core.tools.is_connected

::: neko.core.tools.check_sign

::: neko.core.tools.check_gene_list_format

::: neko.core.tools.mapping_node_identifier

::: neko.core.tools.translate_paths

::: neko.core.tools.consolidate_edges
