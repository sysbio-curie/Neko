# NetworkVisualizer

`NetworkVisualizer` wraps two rendering backends:

- **Graphviz** – static PNG/PDF/SVG diagrams via `render()`
- **yFiles** – interactive widget inside Jupyter via `yfiles_visual()`

## Import

```python
from neko._visual.visualize_network import NetworkVisualizer
```

## Quick example

```python
from neko._visual.visualize_network import NetworkVisualizer

vis = NetworkVisualizer(net, color_by="role")
vis.set_node_colors({"EGFR": "#e74c3c"})
vis.render("my_network", format="svg")
```

---

## Class reference

::: neko._visual.visualize_network.NetworkVisualizer
    options:
      members:
        - __init__
        - set_node_colors
        - set_custom_edge_colors
        - tissue_mapping
        - render
        - yfiles_visual
        - vis_comparison
      show_source: true
      show_root_heading: true
      heading_level: 3
