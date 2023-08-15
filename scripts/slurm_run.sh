#!/bin/bash
#SBATCH -o slurm_run-%A.%a.out
#SBATCH -p RTX3080,GTX980
#SBATCH -n 1
#SBATCH --gres=gpu:1

module load cuda
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"

somd-freenrg -C somd.cfg -p CUDA -m somd.pert -c somd.rst7 -t somd.prm7
