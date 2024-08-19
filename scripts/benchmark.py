import random
import timeit
import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import omnipath as op

from neko.core.network import Network
from neko.inputs import Universe
import itertools

import pywikipathways as pwpw

NETWORKS = dict()

OP_ANNOT_ARGS = {
    'entity_types': 'protein',
    'wide': True,
}

MATRIX = {
    'pathways': 'wikipathway',
    'seeds_percentage': [0.1, 0.2, 0.3],  # Percentage of pathway genes to use as seeds
    'max_len': [2, 3],
    'connection_type': ['complete', 'radial'],
    'network': ['omnipath'], # it does not work if we put multiple resources like omnipath, SIGNOR
}

NUM_ITERATIONS = 5  # Number of iterations for bootstrapping

CUSTOM_PATHWAYS = {
    'EGF/EGFR': 'WP437',  # Large pathway
    'Wnt': 'WP428',  # Medium pathway
    'mTOR': 'WP1471',  # Small pathway
}


#
# Database access
#

def get_pathway_genes(pathway_id):
    return pwpw.get_xref_list(pathway_id, 'H')


def wikipathway_pathways():
    return {
        pathway_name: set(get_pathway_genes(pathway_id))
        for pathway_name, pathway_id in CUSTOM_PATHWAYS.items()
    }


def ensure_msigdb():
    if not isinstance(globals().get('MSIGDB', None), pd.DataFrame):
        globals()['MSIGDB'] = op.requests.Annotations.get(
            resources='MSigDB',
            **OP_ANNOT_ARGS,
        )

    return globals()['MSIGDB'].copy()


def msigdb_pathways(resource: str):
    msigdb = ensure_msigdb()
    msigdb[msigdb.collection == resource]
    msigdb.geneset.cat.remove_unused_categories()
    pathways = dictof_sets(msigdb, 'geneset')

    return pathways


def dictof_sets(df: pd.DataFrame, groupby: str):
    return df.groupby(groupby).genesymbol.apply(set).to_dict()


def other_pathways(resource: str):
    df = op.requests.Annotations.get(resource ** OP_ANNOT_ARGS)

    return dictof_sets(df, 'pathway')


def signalink_pathways():
    return other_pathways('SignaLink_pathway')


def signor_pathways():
    return other_pathways('SIGNOR')


def biocarta_pathways():
    return msigdb_pathways('biocarta_pathways')


def pid_pathways():
    return msigdb_pathways('pid_pathways')


def kegg_pathways():
    return {
        k: v
        for k, v in msigdb_pathways('kegg_pathways').items()
        if (
            'SYNTHESIS' not in k and
            'DEGRADATION' not in k
        )
    }


def ensure_network(resource: str):
    global NETWORKS

    if resource not in NETWORKS:
        NETWORKS[resource] = Universe(resource)

    return NETWORKS[resource]


#
# Benchmark workflow
#

def main():
    result = []

    print("Processing matrix...")

    matrix = MATRIX.copy()
    matrix['iteration'] = range(NUM_ITERATIONS)
    pathways = matrix.pop('pathways')
    pathways = globals()[f'{pathways}_pathways']()
    pathways = {k: v for k, v in pathways.items() if k in CUSTOM_PATHWAYS}
    matrix['pathway'] = pathways.items()

    print("Running analysis...")

    for params in itertools.product(*matrix.values()):
        param = dict(zip(matrix.keys(), params))
        result.append(run(param))

    df = pd.DataFrame(result)
    df.to_csv('network_analysis_results.csv', index=False)

    return df


def run(
    param
):

    network = param['network']
    pathway = param['pathway']
    seed_percentage = param['seeds_percentage']
    max_len = param['max_len']
    connection_type = param['connection_type']
    iteration = param['iteration']

    pw_name, pw_genes = pathway

    print("Processing pathway:", pw_name)

    n_seed = max(3, int(np.round(len(pw_genes) * seed_percentage)))
    seed_genes = random.sample(pw_genes, n_seed)
    param = locals()
    del param['pathway']

    print("Ensuring network...")

    network = ensure_network(network)

    print("Analyzing network...")

    result = analyze_network(
        seed_genes,
        network,
        connection_type,
        max_len,
        pw_genes,
    )

    print("Processing results...")

    result.update(param)
    result['pw_size'] = len(pw_genes)
    result['coverage'] = len(set(result['nodes_found'])) / len(pw_genes)

    return result


def analyze_network(seeds, resources, connection_type, max_len, pw_genes):
    network = Network(seeds, resources=resources.interactions)

    calls = {
        'complete': (
            'network.complete_connection('
            'maxlen=max_len, algorithm="dfs", only_signed=True,'
            'connect_with_bias=False, consensus=False'
            ')'
        ),
        'radial': (
            'network.complete_radially('
            'maxlen=max_len - 1, only_signed=True, consensus=False'
            ')'
        ),
    }

    print("Calling network analysis...")

    exec_time = timeit.timeit(
        stmt=calls[connection_type],
        setup='from __main__ import network, max_len',
        globals=globals(),
        number=1,
    )

    return {
        'num_nodes': len(network.nodes),
        'num_edges': len(network.edges),
        'nodes_found': sorted(pw_genes & set(network.nodes['Genesymbol'])),
        'execution_time': exec_time,
    }

print("Starting analysis...")
final_results = main()

final_results.to_csv("results_benchmark.csv", index=False)

# # Visualization: Multi-variable comparison
# plt.figure(figsize=(20, 12))
# sns.boxplot(x="pathway", y="coverage", hue="seeds_percentage", data=df)
# sns.swarmplot(x="pathway", y="coverage", hue="seeds_percentage", data=df, dodge=True, color=".25", size=3)
#
# # Overlay execution time as color
# for i, pathway in enumerate(df['pathway'].unique()):
#     for j, seed_pct in enumerate(df['seeds_percentage'].unique()):
#         subset = df[(df['pathway'] == pathway) & (df['seeds_percentage'] == seed_pct)]
#         avg_time = subset['execution_time'].mean()
#         plt.scatter(i + (j - 1) * 0.3, subset['coverage'].mean(), s=100, c=[avg_time], cmap='viridis',
#                     vmin=df['execution_time'].min(), vmax=df['execution_time'].max())
#
# plt.colorbar(label='Avg Execution Time (s)')
# plt.title("Pathway Coverage by Seed Percentage with Execution Time")
# plt.savefig("multi_variable_comparison.png")
# plt.close()
#
# # Additional visualizations
#
# # 1. Box plot of coverage by pathway and connection type
# plt.figure(figsize=(15, 10))
# sns.boxplot(x="pathway", y="coverage", hue="connection_type", data=df)
# plt.title("Pathway Coverage by Connection Type")
# plt.savefig("coverage_by_connection_type.png")
# plt.close()
#
# # 2. Scatter plot of number of nodes vs edges, colored by resource type
# plt.figure(figsize=(12, 8))
# sns.scatterplot(x="num_nodes", y="num_edges", hue="resource_type", style="connection_type", size="pathway_size",
#                 data=df)
# plt.title("Number of Nodes vs Edges by Resource Type")
# plt.savefig("nodes_vs_edges_by_resource.png")
# plt.close()
#
# # 3. Line plot of execution time vs max_len
# plt.figure(figsize=(12, 8))
# sns.lineplot(x="max_len", y="execution_time", hue="resource_type", style="connection_type", data=df)
# plt.title("Execution Time vs Max Length")
# plt.savefig("execution_time_vs_max_len.png")
# plt.close()
#
# # 4. Bar plot of average coverage by resource type and pathway
# plt.figure(figsize=(12, 8))
# sns.barplot(x="pathway", y="coverage", hue="resource_type", data=df)
# plt.title("Average Coverage by Resource Type and Pathway")
# plt.savefig("coverage_by_resource_and_pathway.png")
# plt.close()

print("Analysis complete. Results saved to CSV and visualizations generated.")
