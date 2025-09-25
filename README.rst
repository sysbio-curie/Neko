=========================
NeKo: Network Konstructor
=========================

.. figure:: docs/src/neko_logo.png
   :align: right
   :figwidth: 50px
   :alt: NeKo Logo

.. image:: https://github.com/sysbio-curie/Neko/actions/workflows/build.yaml/badge.svg
   :target: https://github.com/sysbio-curie/Neko/actions/workflows/build.yaml
   :alt: Tests

.. image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
   :target: https://sysbio-curie.github.io/Neko/
   :alt: Documentation

Neko is a Python package for extracting, visualizing, converting, and studying interactions from databases into executable activity flow-based models. It's built on top of `Omnipath <https://github.com/saezlab/omnipath>`_, `Pypath <https://github.com/saezlab/pypath>`_, and `Atopo <https://github.com/druglogics/atopo>`_.

Citation
--------

If you use NeKo in your research, please cite our paper:

Ruscone M, Tsirvouli E, Checcoli A, Turei D, Barillot E, et al. (2025) NeKo: A tool for automatic network construction from prior knowledge. *PLOS Computational Biology* 21(9): e1013300. https://doi.org/10.1371/journal.pcbi.1013300

Features
--------

- Network creation and manipulation
- Connection of nodes and subnetworks
- Gene-to-phenotype mapping
- Network visualization
- Interaction database integration

Installation
------------

`NeKo` is still in its alpha version. You can install it from PyPI and also install the necessary external dependencies.

1. **Install `NeKo` from PyPI**:

   Install the main package from PyPI (nekomata, do not confuse with pip install neko or pip install pyneko, those are other packages):

   .. code-block:: bash

       pip install nekomata


Installation from Source
------------------------

For the latest development version, you can still clone the repository and install directly from the source:

.. code-block:: bash

    git clone https://github.com/sysbio-curie/Neko.git
    cd Neko
    pip install .

This will give you the latest version of `NeKo` (not officially released, so be aware there could be some bugs) along with the necessary external dependencies.

Troubleshooting
---------------

If during the installation you encounter problems with the installation of Graphviz, you could be missing basic Graphiz installation on your machine.
You can install it on Linux system with the following command:

.. code-block:: bash

    sudo apt-get install python3-dev graphviz libgraphviz-dev

or on Mac system:

.. code-block:: bash

    brew install python3-dev graphviz libgraphviz-dev

For more details visit: https://graphviz.org/download/

Documentation
-------------

For full documentation, including API reference and detailed tutorials, visit our `GitHub Pages documentation <https://sysbio-curie.github.io/Neko/>`_.

Jupyter Notebooks
-----------------

We provide a comprehensive set of Jupyter notebooks that offer a detailed and user-friendly explanation of the package. These notebooks cover all modules of NeKo and provide a complete overview of how to use the package:


1) Usage
2) Build network using user-defined resources
3) Stepwise connection: a focus on the INE algorithm
4) Connect to upstream components
5) Build network based on kinase-phosphosite interactions
6) Connect to downstream Gene Ontology terms.
7) Map tissue expression
8) Network comparison
9) Re-creating famous pathways from SIGNOR and WIKIPATHWAYS using NeKo


You can find these notebooks in the `notebooks` directory of the repository.

Features comparison with similar tools
--------------------------------------
Below you can find a table displaying the main features of NeKo compared to other similar tools:
`Features Table on GitHub <https://github.com/sysbio-curie/Neko/blob/main/table.md>`_.

Acknowledgements
----------------

This project is a collaborative effort between Institut Curie, NTNU, Saez lab and BSC.

Current contributors: Marco Ruscone, Eirini Tsirvouli, Andrea Checcoli, Dénes Turei, Aasmund Flobak, Emmanuel Barillot, Loredana Martignetti, Julio Saez-Rodriguez and Laurence Calzone.

version 0.9.20
--------------

- Network creation and manipulation: The package allows for the creation of a network of nodes and edges, with various methods for enrichment analysis. This includes adding and removing nodes and edges, loading a network from a SIF (Simple Interaction Format) file, and adding paths to the edge list of the network.
- Database integration: The package provides methods to integrate interactions from databases such as Omnipath, Signor, HURI and others. The user can also integrate personal resource to mine for interactions.
- Database translation: The package provides methods to convert the identifiers of a database storing edges list, into Uniprot.
- Connection of nodes: The package provides several methods to connect nodes in the network. This includes connecting all nodes, connecting a subgroup of nodes, connecting all nodes of a network object, and connecting subcomponents of a network object.
- Connection of genes to phenotype: The package provides a method to connect genes to a phenotype based on provided parameters. This includes retrieving phenotype markers, identifying unique Uniprot genes, and connecting them to the network. There is also an option to compress the network by substituting specified genes with the phenotype name.
