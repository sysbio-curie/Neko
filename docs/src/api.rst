API
===

Import Neko as::

    from neko.core.network import Network

Network
-------

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
