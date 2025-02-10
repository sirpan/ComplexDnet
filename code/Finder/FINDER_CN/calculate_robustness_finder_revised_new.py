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
list_test_data = tuple(new_test[:int(len(new_test))])
test_data_unnumber = pd.DataFrame(list_test_data)
ratio_list =[]
sol_nodes_ratio=[]
sol_nodes_num=[]
sol_nodes_ratio_remain_mean=[]
new_data = pd.DataFrame({'a':test_data_unnumber.iloc[:,0],'b':test_data_unnumber.iloc[:,1],'c':test_data_unnumber.iloc[:,2]})
new_data.to_csv('../data/real/RORGT_finder_renum_temp.txt',sep=' ',index=False,header=None)

#finder
main()

orifinal_data = pd.read_csv('../data/real/sol_temp.csv',header =None).loc[:,0]
list_orifinal_data = list(map(int, orifinal_data))
print(len(list_orifinal_data))
# for ratio in [0.1,0.2]:
for ratio in (0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50):
    sol_nodes_ratio_mean=[]
    sol_nodes_num_mean=[]
    sol_nodes_ratio_remain=[]
    for i in range(2):
        random.shuffle(new_test)
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

        index = []
        index_2 = []
        for i in test_data_unnumber.iloc[:,0]:
            index.append(node_mapping[i])
        for i in test_data_unnumber.iloc[:,1]:
            index_2.append(node_mapping[i])
        new_data = pd.DataFrame({'a':index,'b':index_2,'c':test_data_unnumber.iloc[:,2]})
        new_data.to_csv('../data/real/RORGT_finder_renum_temp.txt',sep=' ',index=False,header=None)

        #finder
        main()

        list_sol = list((pd.read_csv('../data/real/sol_temp.csv',header=None)).loc[:,0])
        sol_new=[]
        for i in list_sol:
            sol_new.append(list(node_mapping.keys())[list(node_mapping.values()).index(i)])
        sol_new = list(map(int, sol_new))
        list_ad_new = list(map(int, list_ad_new))
        # convert str to int
        sol_nodes_num_mean.append(len(set(list_orifinal_data) & set(sol_new)))
        sol_nodes_ratio_mean.append(len(set(list_orifinal_data) & set(sol_new))/len(set(list_orifinal_data)))
        sol_nodes_ratio_remain.append(len(set(list_orifinal_data) & set(sol_new))/len(set(list_orifinal_data)-set(num_nodes_to_remove)))
        sol_nodes_num.append(sol_nodes_num_mean[-1])
        sol_nodes_ratio.append(np.mean(sol_nodes_ratio_mean[-1]))
        sol_nodes_ratio_remain_mean.append(np.mean(sol_nodes_ratio_remain[-1]))
        ratio_list.append(ratio)
# print(sol_nodes_num)
# print(sol_nodes_ratio)
data_final = pd.DataFrame({'ratio':ratio_list,'sol_nodes_num':sol_nodes_num,'sol_nodes_ratio':sol_nodes_ratio,'remain_ratio':sol_nodes_ratio_remain_mean})
data_final.to_csv('robustness_FINDER_block_10.csv',index=False)