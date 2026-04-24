#!/usr/bin/env bash

NPROC=4
BUILD_DIR="/Users/biswajitkhara/projects/mfem-based/mfemprogs/build"

FEXEC="ex0p"
CMD_ARGS="-m ../../data/star.mesh"

COM="mpirun -n ${NPROC} ${BUILD_DIR}/${FEXEC} ${CMD_ARGS}"
echo "${COM}"
${COM} >out.txt 2>&1
