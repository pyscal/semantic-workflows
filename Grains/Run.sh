#!/usr/bin/env bash
# run_all.sh — Atomsk -> LAMMPS -> Slurm, assuming everything is in this folder

set -euo pipefail

#############################
# --- USER CONFIG BLOCK --- #
#############################

# Binaries & files in THIS directory
ATOMSK_BIN="./atomsk"               # must exist in this directory
LAMMPS_BIN= #"/scratch/usr/niihtl18/lammps-29Aug2024/src/lmp_mpi"   # cluster path
POTENTIAL_FILE="Fe_Ack.eam"         # must exist in this directory
LAMMPS_TEMPLATE="in.template.lmp"
SLURM_TEMPLATE="job.template.slurm"

# Geometry / material
LATTYPE="bcc"
ALAT=2.866
ELEM="Fe"
BOX=100                              # Å, cubic box for Voronoi

# Cases
GRAINS=(5 10 15 20)

# Slurm settings
PARTITION="cpu-clx"
WALLTIME="12:00:00"
NNODES=20
TASKS_PER_NODE=96
CPUS_PER_TASK=1
MODULES=(
  "openmpi/gcc/5.0.3"
  "gcc/13.3.0"
)

#########################################
#########################################

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }
BASE_OUT="$(pwd)"

echo "[$(timestamp)] Starting in ${BASE_OUT}"

# 0) Make Fe unit cell once (in cwd)
UNITCELL="Fe_unitcell.xsf"
if [[ ! -f "${UNITCELL}" ]]; then
  echo "[$(timestamp)] Creating unit cell: ${LATTYPE} ${ELEM} (a=${ALAT})"
  "${ATOMSK_BIN}" --create "${LATTYPE}" "${ALAT}" "${ELEM}" "${UNITCELL}"
fi

# 1) Build, write inputs, submit per case
for N in "${GRAINS[@]}"; do
  CASE_DIR="${BASE_OUT}/${N}_grains"
  DATA_FILE="${N}.lmp"

  echo "[$(timestamp)] === N=${N} ==="
  mkdir -p "${CASE_DIR}"
  cp ${POTENTIAL_FILE} ${CASE_DIR}

  # Voronoi params
  VORO="${CASE_DIR}/Fe_voronoi.txt"
  {
    echo "box ${BOX} ${BOX} ${BOX}"
    echo "random ${N}"
  } > "${VORO}"

  # Atomsk polycrystal
  echo "[$(timestamp)] Atomsk -> ${CASE_DIR}/${DATA_FILE}"
  "${ATOMSK_BIN}" --polycrystal "${UNITCELL}" "${VORO}" "${CASE_DIR}/${DATA_FILE}" -wrap
  rm ${CASE_DIR}/${N}_*
	
  # LAMMPS input
  sed -e "s#__DATAFILE__#${DATA_FILE}#g" \
      -e "s#__POTENTIAL__#${POTENTIAL_FILE}#g" \
      "${LAMMPS_TEMPLATE}" > "${CASE_DIR}/in.temp"

  # Slurm job
  MODULE_LOADS_STR=""
  for m in "${MODULES[@]}"; do
    MODULE_LOADS_STR+="module load ${m}\n"
  done

  sed -e "s#__WALLTIME__#${WALLTIME}#g" \
      -e "s#__NNODES__#${NNODES}#g" \
      -e "s#__TPN__#${TASKS_PER_NODE}#g" \
      -e "s#__CPT__#${CPUS_PER_TASK}#g" \
      -e "s#__PARTITION__#${PARTITION}#g" \
      -e "s#__JOBNAME__#${N}#g" \
      -e "s#__LAMMPS_BIN__#${LAMMPS_BIN}#g" \
      -e "s#__MODULE_LOADS__#${MODULE_LOADS_STR%\\n}#g" \
      "${SLURM_TEMPLATE}" > "${CASE_DIR}/job.slurm"

  chmod +x "${CASE_DIR}/job.slurm"

  # Submit
  # echo "[$(timestamp)] sbatch -> ${CASE_DIR}/job.slurm"
  # ( cd "${CASE_DIR}" && sbatch ./job.slurm )
done

echo "[$(timestamp)] All jobs submitted."

