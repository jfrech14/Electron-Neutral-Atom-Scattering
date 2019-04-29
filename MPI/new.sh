#!/bin/bash -l

#SBATCH --job-name= put job name here
#SBATCH --output=%A.txt
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user= put email here
#SBATCH --ntasks=250

enable_lmod
module load icc/19 impi/19 libstdcxx/4

taskset -p -c $$
scontrol show job $SLURM_JOB_ID

mpiexec.hydra ./MPItest
