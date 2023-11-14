#!/bin/bash
#SBATCH -o gmx_run-%A.%a.out
#SBATCH -p main
#SBATCH -n 1
#SBATCH --gres=gpu:1

# module load cuda
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"

srun gmx mdrun -deffnm gromacs -c gromacs_out.gro
