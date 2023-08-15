#!/bin/bash
#SBATCH -o slurm_analysis-%A.%a.out
#SBATCH -p serial
#SBATCH -n 1

srun analyse_freenrg mbar -i lambda*/simfile.dat -p 70 --overlap > dg.txt
