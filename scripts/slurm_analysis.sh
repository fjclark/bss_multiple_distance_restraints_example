#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p serial
#SBATCH -n 1

srun analyse_freenrg mbar -i lambda lambda*/simfile.dat" -p 70 --overlap > dg.txt
