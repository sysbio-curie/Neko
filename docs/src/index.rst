==================
 Neko
==================

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

Neko is a Python package for extracting, visualizing, converting, and studying interactions from databases into executable activity flow-based models. It integrates `OmniPath <https://github.com/saezlab/omnipath>`_ and other interaction resources, uses UniProt tables for identifier translation, and exports networks for tools such as `Atopo <https://github.com/druglogics/atopo>`_.

**Note**: NeKo is distributed as Beta software under the PyPI name
``nekomata``; the Python import package remains ``neko``.

Features
--------

- Network creation and manipulation
- Connection of nodes and subnetworks
- Gene-to-phenotype mapping
- Network visualization and export helpers
- Interaction database integration
- Branching network history with automatic snapshots, HTML/SVG rendering, and pruning controls

Installation
------------

Install the ``nekomata`` distribution from PyPI. Do not confuse it with the
unrelated ``neko`` or ``pyneko`` distributions.

1. **Install `NeKo` from PyPI**:

   .. code-block:: bash

       python -m pip install nekomata

Installation from Source
------------------------

For the latest development version, you can still clone the repository and install directly from the source:

.. code-block:: bash

    git clone https://github.com/sysbio-curie/Neko.git
    cd Neko
    python -m pip install .

This installs the latest development version from the checked-out source.

Troubleshooting
---------------

If Graphviz-related installation or rendering fails, install Graphviz using
your system package manager. On Linux:

.. code-block:: bash

    sudo apt-get install python3-dev graphviz libgraphviz-dev

On macOS:

.. code-block:: bash

    brew install graphviz

For more details visit: https://graphviz.org/download/

Documentation
-------------

For full documentation, including API reference and detailed tutorials, visit our `GitHub Pages documentation <https://sysbio-curie.github.io/Neko/>`_.
Users upgrading an existing workflow should also read the
`NeKo 1.9 migration guide <https://sysbio-curie.github.io/Neko/migration-1.9/>`_.

Jupyter Notebooks
-----------------

We provide a comprehensive set of Jupyter notebooks that offer a detailed and user-friendly explanation of the package. These notebooks cover all modules of NeKo and provide a complete overview of how to use the package:


1) Usage
2) Build network using user-defined resources
3) Stepwise connection: a focus on the INE algorithm
4) Connect to upstream components
5) Build network based on kinase-phosphosite interactions
6) Connect to downstream Gene Ontology terms
7) Map tissue expression
8) Network comparison
9) Re-creating famous pathways from SIGNOR and WIKIPATHWAYS using NeKo
10) Import and complete a network
11) Network history, branching, and visualisation


You can find these notebooks in the `notebooks` directory of the repository.

Features comparison with similar tools
--------------------------------------
Below you can find a table displaying the main features of NeKo compared to other similar tools:
`Features Table on GitHub <https://github.com/sysbio-curie/Neko/blob/main/table.md>`_.

Acknowledgements
----------------

This project is a collaborative effort between Institut Curie, NTNU, Saez lab and BSC.

Current contributors: Marco Ruscone, Eirini Tsirvouli, Andrea Checcoli, Dénes Turei, Aasmund Flobak, Emmanuel Barillot, Loredana Martignetti, Julio Saez-Rodriguez and Laurence Calzone.
