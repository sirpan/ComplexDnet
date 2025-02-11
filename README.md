# NASHnet

This is a framework implementation of NASHnet, as described in our paper:

![demo](https://github.com/sirpan/NASHnet/blob/main/Fig%20.1new.png)

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)


# Overview

Non-alcoholic steatohepatitis (NASH) is a severe form of non-alcoholic fatty liver disease (NAFLD), which is becoming an increasingly prevalent global health issue. Despite ongoing research, effective treatments for NASH remain scarce due to its complex pathogenesis and poorly understood mechanisms. In this study, we developed NASH-predictor, a novel network-based approach that integrates gene expression data and protein-protein interaction networks to identify potential therapeutic targets for NASH. Our method led to the identification of RORγt as a critical target in NASH. We then screened natural compounds and identified Panaxtriol (PXT) from ginseng as a promising candidate. Experimental validation showed that PXT effectively inhibits RORγt and alleviates liver fibrosis in animal models, demonstrating its potential as a therapeutic agent for NASH. These results not only highlight the therapeutic potential of PXT, but also underscore the value of network-based strategies in advancing treatment options for complex diseases like NASH.

# Repo Contents

- [code](./code): source code of NASHnet for the following four cases in the paper.
     - [Finder](./code/Finder): source code for the Finder (FInding key players in complex Networks through DEep Reinforcement learning).
     - [GAPR](./code/GAPR): source code for the GAPR (greedy articulation points removal).
     - [WGCNA](./code/WGCNA): source code for WGCNA (Weighted Gene Co-expression Network Analysis).
     - [gsea](./code/gsea): source code for the ssGSEA(Single-sample gene set enrichment analysis).
     - [SRNA](./code/SRNA): source code for analysis of single-cell RNA-seq data.
- [data](./data): data appled by NASHnet for all code.
- [Results](./Results): Results of paper.


# System Requirements

## Software dependencies and operating systems

### Software dependencies

Users should install the following packages first, which will install in about 5 minutes on a machine with the recommended specs. The versions of software are, specifically:
```
cython==0.29.13 
networkx==2.3 
numpy==1.17.3 
pandas==0.25.2 
scipy==1.3.1 
tensorflow-gpu==1.14.0 
tqdm==4.36.1
```
R dependencies for single-cell genomics analysis, ssGSEA, and WGCNA:
To support single-cell genomics analysis, ssGSEA (single-sample Gene Set Enrichment Analysis), and WGCNA (Weighted Gene Co-expression Network Analysis), install the following R packages:
```
Seurat==4.1.0                 # For single-cell RNA-seq analysis
SingleCellExperiment==1.16.0   # For single-cell data structure
scran==1.22.0                 # For normalization and quality control in single-cell RNA-seq data
scater==1.18.6                # For visualization and analysis of single-cell data
GSVA==1.40.0                  # For Gene Set Variation Analysis (including ssGSEA)
clusterProfiler==4.2.2         # For enrichment analysis and pathway analysis
WGCNA==1.70.3                 # For Weighted Gene Co-expression Network Analysis

```
### Operating systems
The package development version is tested on *Linux and Windows 10* operating systems. The developmental version of the package has been tested on the following systems:

Linux: Ubuntu 18.04  
Windows: 10

The pip package should be compatible with Windows, and Linux operating systems.

Before setting up the FINDER users should have `gcc` version 7.4.0 or higher.

## Hardware Requirements
The model requires a standard computer with enough RAM and GPU to support the operations defined by a user. For minimal performance, this will be a computer with about 4 GB of RAM and 16GB of GPU. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3.3+ GHz/core
GPU: 16+ GB

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores@3.3 GHz) and internet of speed 25 Mbps.






