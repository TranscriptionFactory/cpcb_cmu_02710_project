#!/bin/bash

#SBATCH --job init
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
python py_scripts/initialize.py
