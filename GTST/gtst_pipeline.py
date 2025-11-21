import pandas as pd
import networkx as nx
import GTST
import itertools
import random
import os
import numpy as np

# Ensure consistency
random.seed(123)
np.random.seed(123)


def get_forest(df):
    """
    Convert a forest dataframe into a list of networkx DiGraph objects.
    """
    graphs = []
    for _, row in df.iterrows():
        G = nx.DiGraph()
        n_cases = len(row) + 1  # +1 for root (case 1)

        # Add nodes 1 to N
        G.add_nodes_from(range(1, n_cases + 1))

        # Add edges: row[j] is infector of Case j+2
        # R pipeline: vector index 1 is Case 2.
        for j, infector in enumerate(row):
            case_id = j + 2
            G.add_edge(int(infector), case_id)

        # GTST requires labels; strict string conversion for consistency
        nx.set_node_attributes(G, {n: str(n) for n in G.nodes()}, "label")
        graphs.append(G)
    return graphs


# --- Load Data via Metadata (Exact Matching) ---
meta_path = "Py/data/metadata.csv"
if not os.path.exists(meta_path):
    raise FileNotFoundError("metadata.csv not found. Run R pipeline first.")

metadata = pd.read_csv(meta_path)

# Pre-load forests mapped by tree_id to ensure correct alignment
forest_cache = {}

# We iterate over the metadata to load specific files defined by R
print("Loading forest data...")
for _, row in metadata.iterrows():
    # Reconstruct filename exactly as R wrote it
    fname = f"Py/data/forest_{row['epidemic_size']}_{row['off_R']}_{row['off_k']}_{row['tree_id']}.csv"

    if row["tree_id"] not in forest_cache:
        df = pd.read_csv(fname)
        forest_cache[row["tree_id"]] = get_forest(df)

# Attach forests to metadata for grouping
metadata["forest"] = metadata["tree_id"].map(forest_cache)

# --- Run GTST Analysis ---
sample_sizes = [50, 100]
results = []

print("Starting GTST analysis...")

for es, group in metadata.groupby("epidemic_size"):
    # Sort exactly as R did to ensure ID alignment
    # R: arrange(epidemic_size, off_R, off_k)
    group = group.sort_values(by=["off_R", "off_k"])

    # Convert to list of dicts for easy indexing
    group_records = group.to_dict("records")

    # combinations_with_replacement is equivalent to R's id_A <= id_B
    pairs = list(itertools.combinations_with_replacement(range(len(group_records)), 2))

    for i, j in pairs:
        row_A = group_records[i]
        row_B = group_records[j]

        full_A = row_A["forest"]
        full_B = row_B["forest"]

        for s in sample_sizes:
            # Safety clamp
            curr_s = min(s, len(full_A), len(full_B))

            # Sample without replacement
            forest_A = random.sample(full_A, curr_s)
            forest_B = random.sample(full_B, curr_s)

            # Prepare unique labels for the kernel
            all_nodes = set()
            for G in forest_A + forest_B:
                all_nodes.update(str(n) for n in G.nodes())

            # Initialize and fit GTST
            mmd = GTST.MMD()

            try:
                mmd.fit(
                    G1=forest_A,
                    G2=forest_B,
                    kernel="RW_ARKU_plus",
                    mmd_estimators="MMD_u",
                    r=6,
                    c=0.001,
                    node_label="label",
                    B=999,
                    unique_node_labels=all_nodes,
                    verbose=0,  # Suppress C output if possible
                )
                p_val = mmd.p_values["MMD_u"]
            except Exception as e:
                print(f"Error comparing {row_A['tree_id']} and {row_B['tree_id']}: {e}")
                p_val = None

            results.append(
                {
                    "epidemic_size": es,
                    "tree_id_A": row_A["tree_id"],
                    "off_R_A": row_A["off_R"],
                    "off_k_A": row_A["off_k"],
                    "tree_id_B": row_B["tree_id"],
                    "off_R_B": row_B["off_R"],
                    "off_k_B": row_B["off_k"],
                    "sample_size": curr_s,
                    "p_value": p_val,
                }
            )

results_df = pd.DataFrame(results)
print(results_df.head())
results_df.to_csv("Py/gtst_results.csv", index=False)
