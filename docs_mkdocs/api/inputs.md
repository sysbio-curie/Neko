# Inputs and interaction universes

`Universe` normalizes an interaction source into the DataFrame schema consumed
by `Network`. Built-in adapters cover OmniPath, SIGNOR, HuRI, and
PhosphoSitePlus; a custom DataFrame can be supplied directly.

```python
from neko.inputs import Universe, signor

omnipath_resources = Universe("omnipath")
signor_resources = signor()
custom_resources = Universe(my_interaction_dataframe)
```

A bare `Universe()` is intentionally empty. Use an explicit resource name when
data should be loaded.

## Universe

::: neko.inputs.Universe
    options:
      show_source: true
      show_root_heading: true
      heading_level: 3
      filters:
        - "!^_"

## Adapter functions

::: neko.inputs._universe.network_universe

::: neko.inputs.omnipath

::: neko.inputs.signor

::: neko.inputs.phosphosite

::: neko.inputs.huri
