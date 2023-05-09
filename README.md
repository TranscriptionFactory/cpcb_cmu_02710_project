# jubilant-barnacle
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
