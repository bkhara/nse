#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=4       # number of nodes
#SBATCH --ntasks-per-node=48  # 36 processor core(s) per node
# SBATCH --mem=300G       # maximum memory per node
# SBATCH --constraint=nova21
#SBATCH --job-name="nse"

# Load necessary modules
# module purge
# module load intel
# module load intel-mkl

MPIRUN="mpirun -n ${SLURM_NTASKS}"
BUILD_DIR="${PROJECT_DIR}/build"
FEXEC="nse"
CMD_ARGS="--petscopts ns.petsc --config config.txt"

CWD=`pwd`
cd ${BUILD_DIR} || exit; make ${FEXEC}; cd "${CWD}" || exit;

COM="${MPIRUN} ${BUILD_DIR}/${FEXEC} ${CMD_ARGS}"
echo "${COM}"
${COM}
