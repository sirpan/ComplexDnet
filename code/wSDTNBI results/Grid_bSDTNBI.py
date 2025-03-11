# coding=utf-8
import os
import numpy as np


def makedir(path):
    if os.path.exists(path) is False:
        os.makedirs(path)


def calculate_AUC(V, dataset_1, fingerprint_1, state='process',GPU_command_1='CUDA_VISIBLE_DEVICES=1 '):
    str_V = [str(a) for a in V]
    if state == 'process':
        command =GPU_command_1 + './netinfer -method bnbi -nbi_type SUB DRUG TARGET -nbi_alpha ' + str_V[0] + ' -nbi_beta ' + str_V[1] + ' -nbi_gamma '+str_V[2] + ' -nbi_k 2 -command cv -cv_style nbi -cv_type DRUG TARGET -cv_repeat 1 -cv_fold 10 -cv_seed 12345 -cv_length 5 10 15 20 -training_set DT.tsv+'+fingerprint_1+'.txt -output ./' + dataset_1 + '/CV/' + fingerprint + '/process/CV_OUT_' + ''.join(
            str_V) + fingerprint_1 + '.txt ./' + dataset_1 + '/CV/' + fingerprint_1 + '/process/CV_ROC_' + ''.join(str_V) + fingerprint_1 + '.txt ./' + dataset_1 + '/CV/' + fingerprint_1 + '/process/CV_PR_' + ''.join(str_V) + fingerprint_1 + '.txt'
    elif state == 'end':
        command = GPU_command_1 + './netinfer -method bnbi -nbi_type SUB DRUG Target -nbi_alpha ' + str_V[0] + ' -nbi_beta ' + str_V[1] + ' -nbi_gamma 0' + ' -nbi_k 2 -command cv -cv_style nbi -cv_type DRUG Target -cv_repeat 1 -cv_fold 10 -cv_seed 12345 -cv_length 5 10 15 20 -training_set ../Dataset/Drug_Target_train.txt+../Dataset/' + fingerprint_1 + '_True.txt -output ./' + dataset_1 + '/CV/' + fingerprint + '/result/CV_OUT' + ''.join(
                str_V) + fingerprint_1 + '.txt ./' + dataset_1 + '/CV/' + fingerprint + '/result/CV_ROC' + ''.join(str_V) + fingerprint + '.txt ./' + dataset_1 + '/CV/' + fingerprint + '/result/CV_PR' + ''.join(str_V) + fingerprint + '.txt'
    os.system(command)
    if state == 'process':
        infile = open('./' + dataset_1 + '/CV/' + fingerprint + '/process/CV_OUT_' + ''.join(str_V) + fingerprint_1 + '.txt', 'r')
    else:
        infile = open('./' + dataset_1 + '/CV/' + fingerprint + '/result/CV_OUT_' + ''.join(str_V) + fingerprint_1 + '.txt', 'r')
    AUC = []
    # RE = []
    for index, content in enumerate(infile):
        content = content.rstrip('\n').split('\t')
        if index == 0:
            num_AUC = content.index('AUROC')
            # num_RE = content.index('R_20')
        else:
            AUC.append(float(content[num_AUC]))
            # RE.append(float(content[num_RE]))
    sum_AUC = sum(AUC)
    # sum_RE = sum(RE)
    # V_value = 1 - sum_AUC / len(AUC) + (0.31 - sum_RE / len(RE)) / 6
    # V_value = 0.31 - sum_RE / len(RE)
    V_value = 1 - sum_AUC / len(AUC)
    outfile.write('|'.join(str_V) + '\t' + '%.6f' % V_value + '\n')
    outfile.flush()
    return V_value, str_V


fingerprint_array = ['FCFP4']
# 'PubChem', 'MACCS', 'KR', 'FP4', 'ECFP0', 'ECFP2', , 'FCFP0', 'FCFP2', 'FCFP4'
datasets = ['Grid']
GPU_command = 'CUDA_VISIBLE_DEVICES=1 '
for dataset in datasets:
    for fingerprint in fingerprint_array:
        makedir('./' + dataset + '/CV/' + fingerprint + '/process')
        outfile = open('./' + dataset + '/CV/' + fingerprint + '/AUC_' + fingerprint + '.txt', 'w', encoding='utf-8')
        for x in range(2, 5):
            for y in range(2, 8):
                for z in range(-5,-1):
                    calculate_AUC([round(x * 0.1, 1), round(y * 0.1, 1),round(z*0.1,1)], dataset, fingerprint,GPU_command_1=GPU_command)
