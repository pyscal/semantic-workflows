#!/usr/bin/env bash
#SBATCH -t __WALLTIME__
#SBATCH -N __NNODES__
#SBATCH --tasks-per-node __TPN__
#SBATCH --cpus-per-task=__CPT__
#SBATCH --partition=__PARTITION__
#SBATCH --job-name=__JOBNAME__
#SBATCH --output=slurm-%j.out
#SBATCH --mail-user=__EMAIL__
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
set -x

__MODULE_LOADS__

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export PMI_NO_FORK=1

mpirun -np $SLURM_NTASKS "__LAMMPS_BIN__" -in in.temp
