#!/bin/bash -l

#SBATCH --job-name=Frechem
#SBATCH --output=%A.txt
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jfrec001@odu.edu
#SBATCH --ntasks=125
#SBATCH --exclude=coreV2-22-[001-034]

enable_lmod
module load gcc openmpi

taskset -p -c $$
scontrol show job $SLURM_JOB_ID

mpirun ./MPIScattering
