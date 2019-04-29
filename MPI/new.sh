#!/bin/bash -l

#SBATCH --job-name=Frechem
#SBATCH --output=%A.txt
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jfrec001@odu.edu
#SBATCH --ntasks=5
#SBATCH --constraint=AVX2

enable_lmod
module load icc/19 impi/19 libstdcxx/4

taskset -p -c $$
scontrol show job $SLURM_JOB_ID

mpiexec.hydra ./MPItest
