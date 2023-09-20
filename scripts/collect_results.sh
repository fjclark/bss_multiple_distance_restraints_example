#!/bin/bash
#SBATCH -o slurm_collect_res-%A.%a.out
#SBATCH -p serial
#SBATCH -n 1

# The first argument is the path to the collect_results.py script
srun python $1
