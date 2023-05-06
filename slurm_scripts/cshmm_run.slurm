#!/bin/bash

#SBATCH --job scdiffp
#SBATCH --partition dept_gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mail-user=aar126@pitt.edu
#SBATCH --mail-type=ALL

source ~/.bashrc

job_dir="${USER}_${SLURM_JOB_ID}.dcb.private.net"

mkdir /scr/$job_dir

rsync -a ${SLURM_SUBMIT_DIR} /scr/${job_dir}

trap "echo 'copying files'; rsync -avz * ${SLURM_SUBMIT_DIR}" EXIT

module load anaconda/3-cluster    
module load cuda

eval "$(conda shell.bash hook)"

source activate scdiffenv
python CSHMM-for-time-series-scRNA-Seq-py3/CSHMM_TF_train_release.py --data_file data/tcell_data2.E --structure_file scdiff_runs/tcell_data_filtered_variance.out/tcell_data2.E.h5ad --n_split 20 -ng 1000 --n_iteration 2 --cross_validation 0 --random_seed 0 --model_name test_model_csh
