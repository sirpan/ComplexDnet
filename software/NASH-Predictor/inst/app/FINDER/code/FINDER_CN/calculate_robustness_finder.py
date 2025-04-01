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
    data_test_name = ['2022NADPH_ppi_renum_temp']
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
    os.unlink('../data/real/2022NADPH_ppi_renum_temp.txt')
# get the sol nodes

def main():
    GetSolution(0.01)


test_data = pd.read_table('../data/real/2022NADPH_ppi_renum.txt.',header =None)
new_test = []
for i in range(len(test_data)):
    new_test.append(test_data[0][i].split(' '))
df_new_test = pd.DataFrame(new_test)

orifinal_data = pd.read_table('../data/real/2022NADPH_ppi_renum_CN.txt',header =None)
list_orifinal_data = list(orifinal_data.loc[:,0])
# obtain the original set of sol nodes in FINDER

ratio_list =[]
sol_nodes_ratio=[]
sol_nodes_num=[]



# for ratio in [0.1,0.2]:
for ratio in range(0,99,2):
    ratio=round(ratio*0.01,2)
    random.shuffle(new_test)
    list_test_data = tuple(new_test[:int(len(new_test) * (1.0 - ratio))])
    #shuffle and get the new data to input finder

    # transform the tuple to dataframe
    test_data_unnumber = pd.DataFrame(list_test_data)

    # renumber the nodes from zero to len(data)
    list_ad = list(set(list(test_data_unnumber.iloc[:,0])+(list(test_data_unnumber.iloc[:,1]))))
    list_ad_new = sorted(list_ad)
    list_num = list(range(0,(len(list_ad_new))))
    dict_ad = dict(zip(list_ad_new, list_num))
    index = []
    index_2 = []
    for i in test_data_unnumber.iloc[:,0]:
        index.append(dict_ad[i])
    for i in test_data_unnumber.iloc[:,1]:
        index_2.append(dict_ad[i])
    new_data = pd.DataFrame({'a':index,'b':index_2,'c':test_data_unnumber.iloc[:,2]})
    new_data.to_csv('../data/real/2022NADPH_ppi_renum_temp.txt',sep=' ',index=False,header=None)

    #finder
    main()
    list_sol = list((pd.read_csv('../data/real/sol_temp.csv',header=None)).loc[:,0])
    sol_new=[]
    for i in list_sol:
        sol_new.append(list(dict_ad.keys())[list(dict_ad.values()).index(i)])
    sol_new = list(map(int, sol_new))
    # convert str to int
    sol_nodes_num.append(len(set(list_orifinal_data) & set(sol_new)))
    sol_nodes_ratio.append(len(set(list_orifinal_data) & set(sol_new))/len(list_orifinal_data))
    ratio_list.append(ratio)

# print(sol_nodes_num)
# print(sol_nodes_ratio)
data_final = pd.DataFrame({'ratio':ratio_list,'sol_nodes_num':sol_nodes_num,'sol_nodes_ratio':sol_nodes_ratio})
data_final.to_csv('robustness_2022FINDER.csv',index=False)