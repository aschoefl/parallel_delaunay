#!/bin/bash 
#SBATCH -A edu21.SF2568
#SBATCH -t 00:30:00
#SBATCH --ntasks-per-node=24
module add i-compilers intelmpi
mpirun -np 16 ./delauney
