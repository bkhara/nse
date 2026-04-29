#!/usr/bin/env bash

NPROC=8
BUILD_DIR="${PROJECT_DIR}/build"

FEXEC="nse"
CMD_ARGS="--petscopts ns.petsc --config config.txt"
MPIRUN="mpirun"
#MPIRUN="srun --mpi=pmix"

CWD=`pwd`
cd ${BUILD_DIR} || exit; make ${FEXEC}; cd "${CWD}" || exit;

COM="${MPIRUN} -n ${NPROC} ${BUILD_DIR}/${FEXEC} ${CMD_ARGS}"
echo "${COM}"

${COM} >out.txt 2>&1
