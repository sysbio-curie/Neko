
import random
import time

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import omnipath as op

from neko.core.network import Network
import pywikipathways as pwpw
from neko.inputs import Universe, signor

def get_pathway_genes(pathway_id):
    return pwpw.get_xref_list(pathway_id, 'H')


def analyze_network(seeds, resources, connection_type, max_len):
    start_time = time.time()

    network = Network(seeds, resources=resources)

    if connection_type == "complete":
        network.complete_connection(maxlen=max_len, algorithm="dfs", only_signed=True, connect_with_bias=False,
                                    consensus=False)
    elif connection_type == "radial":
        network.connect_network_radially(max_len=max_len - 1, only_signed=True, consensus=False)

    end_time = time.time()

    return {
        "num_nodes": len(network.nodes),
        "num_edges": len(network.edges),
        "nodes_found": [node for node in all_pathway_genes if node in list(network.nodes["Genesymbol"])],
        "execution_time": end_time - start_time,
    }


# Main script
pathways = {
    "EGF/EGFR": "WP437",  # Large pathway
    "Wnt": "WP428",  # Medium pathway
    "mTOR": "WP1471",  # Small pathway
}

seeds_percentages = [0.1, 0.2, 0.3]  # Percentage of pathway genes to use as seeds
max_lens = [2, 3]
connection_types = ["complete", "radial"]
resource_types = ["Omnipath", "SIGNOR"]
num_iterations = 1  # Number of iterations for bootstrapping

results = []

print("Building resources...")

# Load SIGNOR resources

resources_signor = signor("../neko/_data/signor_db.tsv")  # this function accept only tab separated values
resources_signor.build()
resources_omnipath = Universe()
post_translational = op.interactions.PostTranslational().get()
resources_omnipath.add_resources(post_translational)
resources_omnipath.build()

print("Analyzing networks...")

for pathway_name, pathway_id in pathways.items():
    all_pathway_genes = get_pathway_genes(pathway_id)

    print("Starting analysis for pathway:", pathway_name)

    for seeds_percentage in seeds_percentages: # progressively increase the number of seeds
        seeds_number = max(3, int(len(all_pathway_genes) * seeds_percentage))  # Ensure at least 3 seeds

        print("for seeds_percentage:", seeds_percentage, "seeds_number:", seeds_number)

        for _ in range(num_iterations): # each iteration is a different random seed
            random_seeds = random.sample(all_pathway_genes, seeds_number)

            print("Iteration:", _ + 1)

            for connection_type in connection_types:
                for max_len in max_lens:
                    for resource_type in resource_types:
                        resources = resources_signor.interactions if resource_type == "SIGNOR" else resources_omnipath.interactions

                        network_analysis = analyze_network(random_seeds, resources, connection_type, max_len)

                        results.append({
                            "pathway": pathway_name,
                            "pathway_size": len(all_pathway_genes),
                            "seeds_percentage": seeds_percentage,
                            "seeds_number": seeds_number,
                            "connection_type": connection_type,
                            "max_len": max_len,
                            "resource_type": resource_type,
                            "num_nodes": network_analysis["num_nodes"],
                            "num_edges": network_analysis["num_edges"],
                            "nodes_found": len(network_analysis["nodes_found"]),
                            "coverage": len(network_analysis["nodes_found"]) / len(all_pathway_genes),
                            "execution_time": network_analysis["execution_time"],
                            "initial_seeds": ",".join(random_seeds)
                        })

print("Analysis complete. Processing results...")

# Convert results to DataFrame and save to CSV
df = pd.DataFrame(results)
df.to_csv("network_analysis_results.csv", index=False)

# Visualization: Multi-variable comparison
plt.figure(figsize=(20, 12))
sns.boxplot(x="pathway", y="coverage", hue="seeds_percentage", data=df)
sns.swarmplot(x="pathway", y="coverage", hue="seeds_percentage", data=df, dodge=True, color=".25", size=3)

# Overlay execution time as color
for i, pathway in enumerate(df['pathway'].unique()):
    for j, seed_pct in enumerate(df['seeds_percentage'].unique()):
        subset = df[(df['pathway'] == pathway) & (df['seeds_percentage'] == seed_pct)]
        avg_time = subset['execution_time'].mean()
        plt.scatter(i + (j - 1) * 0.3, subset['coverage'].mean(), s=100, c=[avg_time], cmap='viridis',
                    vmin=df['execution_time'].min(), vmax=df['execution_time'].max())

plt.colorbar(label='Avg Execution Time (s)')
plt.title("Pathway Coverage by Seed Percentage with Execution Time")
plt.savefig("multi_variable_comparison.png")
plt.close()

# Additional visualizations

# 1. Box plot of coverage by pathway and connection type
plt.figure(figsize=(15, 10))
sns.boxplot(x="pathway", y="coverage", hue="connection_type", data=df)
plt.title("Pathway Coverage by Connection Type")
plt.savefig("coverage_by_connection_type.png")
plt.close()

# 2. Scatter plot of number of nodes vs edges, colored by resource type
plt.figure(figsize=(12, 8))
sns.scatterplot(x="num_nodes", y="num_edges", hue="resource_type", style="connection_type", size="pathway_size",
                data=df)
plt.title("Number of Nodes vs Edges by Resource Type")
plt.savefig("nodes_vs_edges_by_resource.png")
plt.close()

# 3. Line plot of execution time vs max_len
plt.figure(figsize=(12, 8))
sns.lineplot(x="max_len", y="execution_time", hue="resource_type", style="connection_type", data=df)
plt.title("Execution Time vs Max Length")
plt.savefig("execution_time_vs_max_len.png")
plt.close()

# 4. Bar plot of average coverage by resource type and pathway
plt.figure(figsize=(12, 8))
sns.barplot(x="pathway", y="coverage", hue="resource_type", data=df)
plt.title("Average Coverage by Resource Type and Pathway")
plt.savefig("coverage_by_resource_and_pathway.png")
plt.close()

print("Analysis complete. Results saved to CSV and visualizations generated.")
