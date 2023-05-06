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
scdiff2 -i newout3.out/tcell_data2.E.h5ad -o scdiff_filterednewout2.out -t /net/dali/home/mscbio/aar126/genomics/scdiff2/tf_dna/Mouse_TF_targets.txt --maxloop 0
#prerun -i tcell_data2.E -o newout3.out
