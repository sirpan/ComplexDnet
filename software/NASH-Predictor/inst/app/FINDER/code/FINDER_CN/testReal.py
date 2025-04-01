#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from FINDER import FINDER
import numpy as np
import time
import pandas as pd
import os

# 启用TensorFlow 1.x兼容模式
# tf.compat.v1.disable_eager_execution()  # 关闭eager execution
# tf = tf.compat.v1
# tf.disable_v2_behavior()

def GetSolution(STEPRATIO):
    dqn = FINDER()
    data_test_path = './'
    data_test_name = ['Finder_PPI']
    model_file_path = './models/'
    model_file_ckpt = 'nrange_30_50_iter_93300.ckpt'
    model_file = os.path.join(model_file_path, model_file_ckpt)  # 使用路径拼接
    save_dir = os.path.join('..', 'results')  # 使用os.path.join生成路径

    # 确保主保存目录存在
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    print('The best model is: %s' % model_file)
    dqn.LoadModel(model_file)
    df = pd.DataFrame(np.arange(1 * len(data_test_name)).reshape((1, len(data_test_name))),
                      index=['time'], columns=data_test_name)
    stepRatio = STEPRATIO

    for j in range(len(data_test_name)):
        print('\nTesting dataset %s' % data_test_name[j])
        data_test = os.path.join(data_test_path, data_test_name[j] + '.txt')
        solution, time = dqn.EvaluateRealData(model_file, data_test, save_dir, stepRatio)
        df.iloc[0, j] = time
        print('Data: %s, time: %.2f' % (data_test_name[j], time))

    # 创建带有StepRatio的子目录
    save_dir_local = os.path.join(save_dir, 'StepRatio_%.4f' % stepRatio)
    if not os.path.exists(save_dir_local):
        os.makedirs(save_dir_local)  # 使用makedirs递归创建目录

    # 保存结果文件
    df.to_csv(os.path.join(save_dir_local, 'sol_time.csv'), encoding='utf-8', index=False)



def EvaluateSolution(STEPRATIO, STRTEGYID):
    #######################################################################################################################
    ##................................................Evaluate Solution.....................................................
    dqn = FINDER()
    data_test_path = './'
    data_test_name = ['Finder_PPI']
    
    # 修正路径生成逻辑，使用os.path拼接
    save_dir = os.path.join('..', 'results', 'StepRatio_%.4f' % STEPRATIO)  # 统一路径格式
    
    # 确保结果目录存在
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    ## begin computing...
    df = pd.DataFrame(np.arange(2 * len(data_test_name)).reshape((2, len(data_test_name))), 
                      index=['solution', 'time'], 
                      columns=data_test_name)
    
    for i in range(len(data_test_name)):
        print('\nEvaluating dataset %s' % data_test_name[i])
        # 使用os.path构建文件路径
        data_test = os.path.join(data_test_path, f"{data_test_name[i]}.txt")
        solution = os.path.join(save_dir, f"{data_test_name[i]}.txt")  # 确保文件名一致
        
        t1 = time.time()
        # strategyID: 0:no insert; 1:count; 2:rank; 3:multiply
        strategyID = STRTEGYID
        score, MaxCCList = dqn.EvaluateSol(data_test, solution, strategyID, reInsertStep=0.001)
        t2 = time.time()
        
        df.iloc[0, i] = score
        df.iloc[1, i] = t2 - t1
        
        # 修正结果文件路径生成
        result_file = os.path.join(save_dir, f'MaxCCList_Strategy_{data_test_name[i]}.txt')
        with open(result_file, 'w') as f_out:
            for val in MaxCCList:
                f_out.write(f'{val:.8f}\n')
        print(f'Data: {data_test_name[i]}, Score: {score:.6f}')
    
    # 保存评估结果到CSV
    result_csv_path = os.path.join(save_dir, 'solution_score.csv')
    df.to_csv(result_csv_path, encoding='utf-8', index=False)

# def main():
    
#     GetSolution(0.01)
#     EvaluateSolution(0.01, 0)


# if __name__=="__main__":
#     main()

