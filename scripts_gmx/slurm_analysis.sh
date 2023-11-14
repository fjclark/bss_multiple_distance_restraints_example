#!/bin/bash
#SBATCH -o gmx_analysis-%A.%a.out
#SBATCH -p main
#SBATCH -n 1

# Was previously using 70 % of the data
srun gmx bar -f lambda*/gromacs.xvg -o -oi -oh > dg.txt