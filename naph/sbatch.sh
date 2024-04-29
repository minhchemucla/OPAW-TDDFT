#!/bin/bash
#SBATCH --job-name="0.5"
#SBATCH --partition=shared
#SBATCH --constraint="lustre"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -c 1
#SBATCH --account=cla301
#SBATCH --export=ALL
#SBATCH -t 24:00:00

export SLURM_EXPORT_ENV=ALL
#export OMP_NUM_THREADS=4

module load cpu/0.17.3b gcc/10.2.0/npcyll4 intel-mkl/2020.4.304/ghfk3mu openmpi/4.1.1/ygduf2r

mpirun -n 32 ./opaw_tddft.x > log 2> err 
