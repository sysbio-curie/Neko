API Reference
=============

This page contains the API reference for the Neko package.

Network Class
-------------

The main class in NeKo is the ``Network`` class. You can import it as follows:

.. code-block:: python

    from neko.core.network import Network

Below is a summary of the methods available in the ``Network`` class:

.. autoclass:: neko.core.network.Network
   :members:
   :undoc-members:
   :show-inheritance:

Method Details
--------------

.. autosummary::
   :toctree: _autosummary

   neko.core.network.Network.load_network_from_sif
   neko.core.network.Network.add_edge
   neko.core.network.Network.remove_node
   neko.core.network.Network.add_cascade_to_edge_list
   neko.core.network.Network.add_node
   neko.core.network.Network.add_paths_to_edge_list
   neko.core.network.Network.complete_connection
   neko.core.network.Network.connect_component
   neko.core.network.Network.connect_genes_to_phenotype
   neko.core.network.Network.connect_nodes
   neko.core.network.Network.connect_subgroup
   neko.core.network.Network.connect_to_upstream_nodes
   neko.core.network.Network.convert_edgelist_into_genesymbol
   neko.core.network.Network.filter_unsigned_paths
   neko.core.network.Network.is_connected

NetworkVisualizer Class Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: _autosummary

   neko._visual.visualize_network.NetworkVisualizer.set_custom_edge_colors
   neko._visual.visualize_network.NetworkVisualizer.set_node_colors
   neko._visual.visualize_network.NetworkVisualizer.add_edges_to_graph
   neko._visual.visualize_network.NetworkVisualizer.add_nodes_to_graph
   neko._visual.visualize_network.NetworkVisualizer.tissue_mapping
   neko._visual.visualize_network.NetworkVisualizer.render
   neko._visual.visualize_network.NetworkVisualizer.yfiles_visual
   neko._visual.visualize_network.NetworkVisualizer.vis_comparison

Other Modules
-------------

.. Add other modules and classes here as your project grows
