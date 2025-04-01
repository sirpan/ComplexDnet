# Written by Zengrui Wu, LMMD, ECUST, 2017-10-25
# Written by Ze Wang, LMMD, ECUST, 2022-11-2

import pandas as pd
import numpy as np
import random, sys, os
from FINDER import FINDER

dqn = FINDER()
model_file_path = './models/'
model_file_ckpt = 'nrange_30_50_iter_93300.ckpt'
model_file = model_file_path + model_file_ckpt
dqn.LoadModel(model_file)
def GetSolution(STEPRATIO):
    ######################################################################################################################
    ##................................................Get Solution (model).....................................................
    data_test_path = '../data/real/'
    data_test_name = ['RORGT_finder_renum_temp']
    ## save_dir
    save_dir = '../results/FINDER_CN/test'
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    #################################### modify to choose which stepRatio to get the solution
    stepRatio = STEPRATIO
    for j in range(len(data_test_name)):
        print ('\nTesting dataset %s'%data_test_name[j])
        data_test = data_test_path + data_test_name[j] + '.txt'
        solution, time = dqn.EvaluateRealData(model_file, data_test, save_dir, stepRatio)
    pd.DataFrame(solution).to_csv('../data/real/sol_temp.csv',index=False,header=None)
    os.unlink('../data/real/RORGT_finder_renum_temp.txt')
# get the sol nodes

def main():
    GetSolution(0.01)


test_data = pd.read_table('../data/real/block_new_ppi_renum.txt',header =None)
new_test = []
for i in range(len(test_data)):
    new_test.append(test_data[0][i].split(' '))
df_new_test = pd.DataFrame(new_test)

orifinal_data = pd.read_table('../data/real/2022RORGT_block_ppi_CN.txt',header =None)
list_orifinal_data = list(orifinal_data.loc[:,0])
# obtain the original set of sol nodes in FINDER

ratio_list =[]
sol_nodes_ratio=[]
sol_nodes_num=[]
sol_nodes_ratio_remain_mean=[]


# for ratio in [0.1,0.2]:
for ratio in np.arange(0.01, 0.51, 0.01):  # ratios from 1% to 50%
    sol_nodes_ratio_mean = []
    sol_nodes_num_mean = []
    sol_nodes_ratio_remain = []

    # Repeat the experiment multiple times for robustness
    for _ in range(2):  # Adjust the repeat count as needed for robustness
        random.shuffle(new_test)
        
        # Determine the number of nodes to remove
        num_nodes_to_remove = int(len(new_test) * ratio)
        num_nodes_to_remove = max(num_nodes_to_remove, 1)  # Ensure at least one node is removed

        # Randomly select nodes to remove
        nodes = sorted(set([node for pair in new_test for node in pair]))  # Extract unique nodes
        nodes_to_remove = random.sample(nodes, num_nodes_to_remove)

        # Create the reduced dataset by removing the selected nodes
        remaining_data = [pair for pair in new_test if pair[0] not in nodes_to_remove and pair[1] not in nodes_to_remove]
        test_data_unnumber = pd.DataFrame(remaining_data)

        # Renumber the remaining nodes
        remaining_nodes = sorted(set([node for pair in remaining_data for node in pair]))
        node_mapping = {node: idx for idx, node in enumerate(remaining_nodes)}
        print(test_data_unnumber.shape[0])
        index = []
        index_2 = []
        for i in test_data_unnumber.iloc[:, 0]:
            index.append(node_mapping[i])
        for i in test_data_unnumber.iloc[:, 1]:
            index_2.append(node_mapping[i])

        new_data = pd.DataFrame({'a': index, 'b': index_2, 'c': test_data_unnumber.iloc[:, 2]})
        new_data.to_csv('../data/real/RORGT_finder_renum_temp.txt', sep=' ', index=False, header=None)

        # Call FINDER to evaluate the modified data
        main()

        # Read the solution nodes produced by FINDER
        list_sol = list(pd.read_csv('../data/real/sol_temp.csv', header=None).iloc[:, 0])
        sol_new = [list(node_mapping.keys())[list(node_mapping.values()).index(i)] for i in list_sol]
        sol_new = list(map(int, sol_new))
        #nodes_to_remove = list(map(int, nodes_to_remove))

        # Calculate metrics based on original and remaining nodes
        intersection = set(list_orifinal_data) & set(sol_new)

        remaining_intersection = set(list_orifinal_data) & set(remaining_nodes)

        # Append the calculated values to the corresponding lists
        sol_nodes_num_mean.append(len(intersection))
        sol_nodes_ratio_mean.append(len(intersection) / len(set(list_orifinal_data)))
        sol_nodes_ratio_remain.append(len(intersection) / len(remaining_intersection) if len(remaining_intersection) > 0 else 0)
        sol_nodes_num.append(sol_nodes_num_mean[-1])
        sol_nodes_ratio.append(np.mean(sol_nodes_ratio_mean[-1]))
        sol_nodes_ratio_remain_mean.append(np.mean(sol_nodes_ratio_remain[-1]))
        ratio_list.append(ratio)
# print(sol_nodes_num)
# print(sol_nodes_ratio)
data_final = pd.DataFrame({'ratio':ratio_list,'sol_nodes_num':sol_nodes_num,'sol_nodes_ratio':sol_nodes_ratio,'remain_ratio':sol_nodes_ratio_remain_mean})
data_final.to_csv('robustness_FINDER_block_10.csv',index=False)