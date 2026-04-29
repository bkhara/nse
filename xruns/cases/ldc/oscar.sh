#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -J dyb
#SBATCH -o out.%j    # Name of stdout output file
#SBATCH -e out.%j    # Name of stderr error file

#SBATCH -p batch     # Queue (partition) name
#SBATCH -t 6:00:00  # Run time (hh:mm:ss)
#SBATCH --nodes=2         # Total #of nodes
#SBATCH --tasks-per-node=32        # Total processes per node

module load hpcx-mpi/2.25.1s-le4f

OSCAR_MPIRUN="srun --mpi=pmix"
BUILD_DIR="${PROJECT_DIR}/build"
FEXEC="frac"
CMD_ARGS="--petscopts frac.petsc --config config.txt"
#CMD_ARGS="--petscopts frac.petsc.vi --config config.txt"
#CMD_ARGS="--petscopts frac.petsc.bddc --config config.txt"

COM="${OSCAR_MPIRUN} ${BUILD_DIR}/${FEXEC} ${CMD_ARGS}"
echo "${COM}"
${COM}
