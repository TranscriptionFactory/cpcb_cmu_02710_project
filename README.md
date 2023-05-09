# jubilant-barnacle
## Getting Data
use gdown https://drive.google.com/drive/folders/1XiVwOaTIxERE5Kf2T2lAfuW8sLL3IS8u?usp=sharing to download files from paper


## Analysis Notebooks
Analysis ipynbs show results from CSHMM training for various different runs 

- scRNA-seq data only: (cv = cross validation, kc6 = starting with 6 initial clusters in scdiff, split6 = 6 initial splits during CSHMM training
- scRNA-seq/scATAC-seq: (integrated)


## Other Notebooks and scripts
- monocle3_analysis.R: processing steps for analyzing data with monocle
- format_data.ipynb: formatting raw data & renaming clusters into pseudotime 
- get_paper_cluster_markers.ipynb: getting top (by p-value) markers for each cluster 
- run_cshmm.ipynb: testing running CSHMM training (before moving to doing training on the cluster)
- output_analysis.ipynb: exploratory data analysis for CSHMM training outputs


## Setup: 
### conda_envs folder
- conda_envs/env_setup.ipynb: instructions for installing dependencies for cshmm & scdiff2 (must be done in separate conda envs because they use different versions of python)
- final_cshmm_env.yml: conda env spec for running cshmm
- scdiff2_env: conda env spec for running scdiff2

### data processing
- format_data.ipynb: processing raw scRNA-seq data for CSHMM and scdiff
- formated_integrated_data.ipynb: processing integrated scRNA-seq/scATAC-seq for CSHMM and scdiff
