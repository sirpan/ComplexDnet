# Written by Zengrui Wu, LMMD, ECUST, 2017-10-25
# Written by Ze Wang, LMMD, ECUST, 2022-11-2
import pandas as pd
import os
import numpy as np
import random
import networkx as nx
from FINDER import FINDER

# Initialize the FINDER model
dqn = FINDER()
model_file_path = './models/'
model_file_ckpt = 'nrange_30_50_iter_93300.ckpt'
model_file = model_file_path + model_file_ckpt
dqn.LoadModel(model_file)

# Define the function to evaluate data using FINDER
def GetSolution(STEPRATIO):
    data_test_path = '../data/real/'
    data_test_name = ['RORGT_finder_renum_temp']
    save_dir = '../results/FINDER_CN/test'
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    stepRatio = STEPRATIO
    for j in range(len(data_test_name)):
        data_test = data_test_path + data_test_name[j] + '.txt'
        solution, time = dqn.EvaluateRealData(model_file, data_test, save_dir, stepRatio)
    pd.DataFrame(solution).to_csv('../data/real/sol_temp.csv', index=False, header=None)
    os.unlink('../data/real/RORGT_finder_renum_temp.txt')

# Function to run the FINDER model
def main():
    GetSolution(0.01)

# Load the graph data
df = pd.read_csv('../data/real/block_new_ppi_renum.txt', sep=' ',header=None)
G_original = nx.from_pandas_edgelist(df, 0, 1)
original_APs = set(nx.articulation_points(G_original))
original_nodes = set(G_original.nodes())

orifinal_data = pd.read_table('../data/real/2022RORGT_block_ppi_CN.txt',header =None)
list_orifinal_data = list(orifinal_data.loc[:,0])
# Prepare lists to store results
percentages = np.arange(1, 10)  # 1% to 50%
sol_nodes_num = []
sol_nodes_ratio = []
sol_nodes_ratio_remain_mean = []
ratio_list = []

# Iterate over each removal percentage
for pct in percentages:
    # Calculate the number of nodes to remove
    num_nodes_to_remove = int(len(G_original) * pct / 100)
    num_nodes_to_remove = min(num_nodes_to_remove, len(G_original))  # Ensure it doesn't exceed the number of nodes

    # Randomly select nodes to remove
    nodes_to_remove = random.sample(original_nodes, num_nodes_to_remove)

    # Create a copy of the graph and remove the selected nodes
    G_removed = G_original.copy()
    G_removed.remove_nodes_from(nodes_to_remove)

    # Extract remaining edges after node removal
    remaining_data = [(u, v) for u, v in G_removed.edges()]

    # Prepare the data for FINDER
    test_data_unnumber = pd.DataFrame(remaining_data, columns=['a', 'b'])

    remaining_nodes = sorted(set(test_data_unnumber['a'])|(set(test_data_unnumber['b'])))
    node_mapping = {node: idx for idx, node in enumerate(remaining_nodes)}

    # Renumber the nodes and save the renumbered data for FINDER input
    index = [node_mapping[u] for u in test_data_unnumber['a']]
    index_2 = [node_mapping[v] for v in test_data_unnumber['b']]
    new_data = pd.DataFrame({'a': index, 'b': index_2, 'c': "{}"})
    new_data.to_csv('../data/real/RORGT_finder_renum_temp.txt', sep=' ', index=False, header=None)

    # Evaluate the network with FINDER
    main()

    # Read the solution nodes from FINDER
    list_sol = pd.read_csv('../data/real/sol_temp.csv', header=None).iloc[:, 0]
    
    
    sol_new = [list(node_mapping.keys())[list(node_mapping.values()).index(i)] for i in list_sol]
    # Calculate metrics based on original and remaining nodes
    intersection = set(list_orifinal_data) & set(sol_new)
    remaining_intersection = set(list_orifinal_data) & set(remaining_nodes)

    # Append results to the corresponding lists
    sol_nodes_num.append(len(intersection))
    sol_nodes_ratio.append(len(intersection) / len(set(list_orifinal_data)))
    sol_nodes_ratio_remain_mean.append(len(intersection) / len(remaining_intersection) if len(remaining_intersection) > 0 else 0)
    ratio_list.append(pct)

# Save the results to a DataFrame and export to CSV
results_df = pd.DataFrame({
    'Ratio': ratio_list,
    'Solution Nodes Count': sol_nodes_num,
    'Solution Nodes Ratio': sol_nodes_ratio,
    'Remain Ratio': sol_nodes_ratio_remain_mean
})
results_df.to_csv('robustness_FINDER_block_10.csv', index=False)
