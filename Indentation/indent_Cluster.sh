#!/bin/bash
#SBATCH -t 12:00:00 # walltime
#SBATCH -N 60 # number of nodes
#SBATCH --tasks-per-node 96
#SBATCH -p standard96
#SBATCH --mail-user=hoang-thien.luu@tu-clausthal.de # email address
#SBATCH --mail-type=BEGIN,END,FAIL

set -x
module load impi/2019.5
module load gcc/9.3.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export THISDIR=$(pwd)
export PMI_NO_FORK=1
export LAMMPS_POTENTIALS=${THISDIR}

set -x
export THISDIR=$(pwd)

Vel=('V26' 'V28' 'V30' 'V32' 'V34' 'V36' 'V38' 'V40' 'V42' 'V44' 'V46' 'V48' 'V50' 'V75' 'V100' 'V125' 'V150' 'V175' 'V200')
Vel2=('0.26' '0.28' '0.3' '0.32' '0.34' '0.36' '0.38' '0.4' '0.42' '0.44' '0.46' '0.48' '0.5' '1' '0.75' '1.25' '1.5' '1.75' '2')
Dum=('27692' '25714' '24000' '22500' '21176' '20000' '18947' '18000' '17142' '16262' '15652' '15000' '14400' '7200' '9600' '5760' '4800' '4110' '3600')
Pri=('46' '42' '40' '37' '35' '34' '32' '30' '28' '28' '26' '25' '24' '12' '16' '10' '8' '7' '6')

###################### ${THISDIR} ####################
####################### Level 1 #########################
for (( i = 0; i < ${#Vel[@]};i++ ));  do
  #### Variable 
  Velocity=${Vel[$i]}
  Velocity2=${Vel2[$i]}
  Print=${Pri[$i]}
  Dump=${Dum[$i]}
  mkdir ${Vel[$i]}
  cd ${Vel[$i]}
  cp ${THISDIR}/in.indent ${THISDIR}/${Vel[$i]}/
  mpirun -np $SLURM_NTASKS /scratch/usr/nicngu18/lammps-3Mar20/src/lmp_mpi -var Pri ${Pri[$i]} -var Dum ${Dum[$i]} -var Vel2 ${Vel2[$i]} -in in.indent &
  wait
  ###################### End Level 1 ####################
  cd ..
###################### ${THISDIR} ####################
done
exit
