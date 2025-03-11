#!/bin/bash
#SBATCH --gpus=1
./netinfer -method wnbi -nbi_type SUB COMPOUND+DRUG TARGET -nbi_k 2 -nbi_alpha 0.4 -nbi_beta 0.2 -nbi_gamma -0.5 ©\nbi_delta 20 ©\nbi_epsilon 4 -command predict -predict_type COMPOUND+DRUG TARGET -predict_num 50 -training_set DT.tsv+FCFP4.txt+FCFP4_nature.txt -output result_nature.txt


