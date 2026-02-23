# Tutorials

The tutorials below are Jupyter notebooks that walk you through every major NeKo feature, from building your first network to comparing two signalling models.

| # | Notebook | Topics covered |
|---|---|---|
| 1 | [Network Building](1_network_building.ipynb) | Create a `Network`, load resources, `connect_nodes` |
| 2 | [Add Resources](2_add_resources.ipynb) | Extend the interaction universe with custom databases |
| 3 | [Stepwise Connection](3_stepwise_connection.ipynb) | Fine-grained connection strategies |
| 4 | [Connect Upstream](4_Connect_upstream.ipynb) | Link upstream regulators to a seed node set |
| 5 | [Phosphosite Network](5_build_phosphosite_network.ipynb) | Build a kinase–substrate network from PhosphoSitePlus |
| 6 | [Ontology](6_ontology.ipynb) | Map nodes to Gene Ontology terms and phenotypes |
| 7 | [Tissue Mapping](7_tissue_mapping.ipynb) | Filter interactions by tissue expression |
| 8 | [Compare Networks](8_Compare_networks.ipynb) | Diff two networks and highlight differences |
| 9 | [Recreating Famous Pathways](9-Recreating_famous_pathways.ipynb) | Reproduce well-known signalling diagrams |
| 10 | [Import & Complete a Network](10_Import_and_complete_a_network.ipynb) | Import from SIF / GraphML and fill gaps |
| 11 | [Network History](11_network_history.ipynb) | Snapshot, branch, and diff network states |

---

## Running the notebooks locally

```bash
git clone https://github.com/sysbio-curie/Neko.git
cd Neko
pip install nekomata
jupyter lab notebooks/
```

!!! note
    Notebooks that download data (e.g., from OmniPath) require internet access and might take a minute on first run while the database is cached locally.
