==================
Neko
==================

.. image:: https://github.com/sysbio-curie/Neko/actions/workflows/build.yaml/badge.svg
   :target: https://github.com/sysbio-curie/Neko/actions/workflows/build.yaml
   :alt: Tests

.. image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
   :target: https://sysbio-curie.github.io/Neko/
   :alt: Documentation

Neko: Network Konstructor
-------------------------------------------------------------------

Neko is a Python package for extracting, visualizing, converting, and studying interactions from databases into executable activity flow-based models. It's built on top of `Omnipath <https://github.com/saezlab/omnipath>`_, `Pypath <https://github.com/saezlab/pypath>`_, and `Atopo <https://github.com/druglogics/atopo>`_.

**Note**: Neko is currently in development and approaching its final stages. It is not yet available on PyPI.

Features
--------

- Network creation and manipulation
- Connection of nodes and subnetworks
- Gene-to-phenotype mapping
- Network visualization
- Interaction database integration

Installation
------------

As Neko is still in development, you can install it directly from the GitHub repository:

.. code-block:: bash

    git clone https://github.com/sysbio-curie/Neko.git
    cd Neko
    pip install -e .

Documentation
-------------

For full documentation, including API reference and detailed tutorials, visit our `GitHub Pages documentation <https://sysbio-curie.github.io/Neko/>`_.

Jupyter Notebooks
-----------------

We provide a comprehensive set of Jupyter notebooks that offer a detailed and user-friendly explanation of the package. These notebooks cover all modules of NeKo and provide a complete overview of how to use the package:

1. Network Building
2. Adding Resources
3. Building Phosphosite Networks
4. Connecting Upstream
5. Ontology
6. Tissue Mapping

You can find these notebooks in the `notebooks` directory of the repository.

Acknowledgements
----------------

This project is a collaborative effort with Dénes Turei and Asmund Flobak.

Current contributors: Marco Ruscone, Eirini Tsirvouli, Andrea Checcoli, Dénes Turei.

version 0.9.1
--------------

- Network creation and manipulation: The package allows for the creation of a network of nodes and edges, with various methods for enrichment analysis. This includes adding and removing nodes and edges, loading a network from a SIF (Simple Interaction Format) file, and adding paths to the edge list of the network.
- Connection of nodes: The package provides several methods to connect nodes in the network. This includes connecting all nodes, connecting a subgroup of nodes, connecting all nodes of a network object, and connecting subcomponents of a network object.
- Connection of genes to phenotype: The package provides a method to connect genes to a phenotype based on provided parameters. This includes retrieving phenotype markers, identifying unique Uniprot genes, and connecting them to the network. There is also an option to compress the network by substituting specified genes with the phenotype name.
