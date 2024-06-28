==================
Neko
==================

.. image:: https://github.com/sysbio-curie/Neko/actions/workflows/test.yml/badge.svg
   :target: https://github.com/sysbio-curie/Neko/actions/workflows/test.yml
   :alt: Tests

.. image:: https://img.shields.io/readthedocs/omniflow_project
   :target: https://neko.readthedocs.io
   :alt: Documentation

About Neko
----------

Neko is a package to extract, visualize, convert and study interactions from databases into executable activity flow based models.
It is based on `Omnipath <https://github.com/saezlab/omnipath>`_, `Pypath <https://github.com/saezlab/pypath>`_ and `Atopo <https://github.com/druglogics/atopo>`_.

This is a work-in-progress package made in collaboration with Dénes Turei and Asmund Flobak.

Current contributors: Marco Ruscone, Eirini Tsirvouli, Andrea Checcoli, Dénes Turei.

Installation
------------

All the dependencies are listed in the pyproject.toml file. To install the package, you can use the following command:

.. code-block:: bash

    pip install .

Current Features
----------------

version 0.9.0
^^^^^^^^^^^^^

- **Network creation and manipulation**:
    - Create networks of nodes and edges
    - Add and remove nodes and edges
    - Load networks from SIF (Simple Interaction Format) files
    - Add paths to the edge list of the network
    - Various methods for enrichment analysis

- **Connection of nodes**:
    - Connect all nodes
    - Connect a subgroup of nodes
    - Connect all nodes of a network object
    - Connect subcomponents of a network object

- **Connection of genes to phenotype**:
    - Connect genes to a phenotype based on provided parameters
    - Retrieve phenotype markers
    - Identify unique Uniprot genes and connect them to the network
    - Option to compress the network by substituting specified genes with the phenotype name

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

   contents
   api
