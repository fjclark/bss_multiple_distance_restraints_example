#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p GTX980,RTX3080
#SBATCH -n 1
#SBATCH --gres=gpu:1

module load cuda
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"

srun somd-freenrg -C somd.cfg -p CUDA -m somd.pert -c somd.rst7 -t somd.prm7
