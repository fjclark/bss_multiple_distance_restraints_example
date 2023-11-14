#!/bin/bash
#SBATCH -o gmx_collect_res-%A.%a.out
#SBATCH -p main
#SBATCH -n 1

# The first argument is the path to the collect_results.py script
srun python $1
