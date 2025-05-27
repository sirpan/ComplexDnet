# ComplexDnet

This is a framework implementation of ComplexDnet, as described in our paper:

![demo](https://github.com/sirpan/ComplexDnet/blob/main/Fig%20.1new.png)

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [User guide](#User_guide)
# Overview

The intricate pathogenesis of complex diseases poses formidable challenges for target and drug discovery, particularly in disorders like non-alcoholic steatohepatitis (NASH) where cross-talking multi-system induced insufficient efficacy of treatment.  Here, we present ComplexDnet—a first-in-class AI-driven framework integrating transcriptomic data, protein-protein (PPI) network to prioritize targets and decode mechanisms. Unlike conventional PPI-centric approaches, ComplexDnet establishes a transcriptome-driven framework that identifies different types of complex disease novel targets, achieving average 77.63% recall in eight types of cancer target identification, surpassing traditional differential gene (DEG) analysis by ~20%. Then we apply ComplexDnet in NASH and revealing RORγt as a central regulator of NASH-associated inflammatory and fibrotic cascades. Leveraging this finding, we conducted a network-based VS on RORγt, and 28 compounds were selected for experimental validation. In a high-throughout RORγt-GAL4 luciferase reporter screen, panaxatriol (PXT) was identified as a novel, direct RORγt inhibitor (IC50=0.01 µM). The 2.8 Å resolution X-ray structure of the RORγt-PXT complex confirmed PXT acts as an inverse agonist. Given RORγt's role in NASH-driven liver fibrosis, we evaluated PXT in mouse liver fibrosis models, where it significantly attenuated collagen deposition and inflammatory responses. Critically, we packaged this pipeline into an open-source software (https://github.com/sirpan/ComplexDnet ), achieving higher hit rates and speed than random screening in target identification. ComplexDnet provides a AI-driven framework for NASH therapeutic strategy discovery, with a transferable framework for multi-disease treatment optimization.

# Repo Contents

- [code](./code): source code of ComplexDnet for the following four cases in the paper.
     - [Finder](./code/Finder): source code for the Finder (FInding key players in complex Networks through DEep Reinforcement learning).
     - [GAPR](./code/GAPR): source code for the GAPR (greedy articulation points removal).
     - [WGCNA](./code/WGCNA): source code for WGCNA (Weighted Gene Co-expression Network Analysis).
     - [gsea](./code/gsea): source code for the ssGSEA(Single-sample gene set enrichment analysis).
     - [SRNA](./code/SRNA): source code for analysis of single-cell RNA-seq data.
     - [wSDTNBI results](./code/wSDTNBI%20results): source code and results for prediction of DTI using wSDTNBI.
     - [figure code](./code/figure%20code): source code and results for figure.
     - [NASH-Score](./code/NASH-Score): source code and results for NASH-Score.
- [data](./data): data appled by ComplexDnet for all code.
- [Results](./Results): Results of paper.
- [software](./software): software of NASH-Predictor.

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

# User_guide
Before using the NASH-Predictor software, please read the Directions for use.pdf in the root directory carefully




