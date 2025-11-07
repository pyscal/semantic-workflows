#!/bin/bash
#SBATCH --job-name=fe_bcc_hcp_transition
#SBATCH --time=23:59:00
#SBATCH --partition=cpu
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=3GB
#SBATCH --hint=nomultithread
#SBATCH --chdir=/home/menonsqr/repos/FeC-application/phase_transition
#SBATCH --account=drautrmy_0004
#SBATCH -N 1
ulimit -l unlimited
source activate fec2
export UCX_MEMTYPE_CACHE=n
export UCX_ERROR_SIGNALS=""
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp

python bcc_hcp.py
