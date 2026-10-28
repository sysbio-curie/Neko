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

## Custom interaction tables

A custom DataFrame needs `source` and `target` columns. Give interaction signs
either as boolean `is_stimulation` / `is_inhibition` columns or as one lowercase
`effect` column. Recognised `effect` values (case-insensitive) are:

| Meaning | Values |
| --- | --- |
| Stimulation | `stimulation`, `activation`, `activate(s)`, `stimulate(s)`, `positive`, `+`, `True`, `1` |
| Inhibition | `inhibition`, `inhibit(s)`, `repression`, `repress(es)`, `negative`, `-`, `False`, `-1` |
| Both | `bimodal`, `both` |
| Undefined | `undefined`, `unknown`, or an empty value |
| Complex formation | `form complex`, `complex` |

Any other value raises `ValueError` and names the offending values; NeKo does
not guess a sign. Boolean flag columns (`is_directed`, `is_stimulation`,
`is_inhibition`, `form_complex`) accept `True`/`False`, `yes`/`no`, `t`/`f`,
`y`/`n` and `1`/`0`/`-1` in any case, and raise `ValueError` for other text.

```python
import pandas as pd
from neko.core.network import Network

table = pd.DataFrame({
    "source": ["EGFR", "PTEN"],
    "target": ["PIK3CA", "AKT1"],
    "effect": ["stimulation", "inhibition"],
})
net = Network(["EGFR", "AKT1"], resources=Universe(table).interactions)
```

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
