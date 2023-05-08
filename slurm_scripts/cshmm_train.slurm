#!/bin/bash

#SBATCH --job cstr
#SBATCH --partition dept_gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mail-user=aar126@pitt.edu


source ~/.bashrc

job_dir="${USER}_${SLURM_JOB_ID}.dcb.private.net"

mkdir /scr/$job_dir

rsync -a ${SLURM_SUBMIT_DIR} /scr/${job_dir}

trap "echo 'copying files'; rsync -avz * ${SLURM_SUBMIT_DIR}" EXIT

module load anaconda/3-cluster    
module load cuda

eval "$(conda shell.bash hook)"

source activate cshmm
python CSHMM-for-time-series-scRNA-Seq/CSHMM_train.py --data_file data/integrated_data.E --structure_file data/integrated_data.Einit_cluster.txt --n_split 100 -ng 10000 --n_iteration 20 --cross_validation 0 --random_seed 0 --model_name integrated_data