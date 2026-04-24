force_rebuild=false

BUILD_TYPE="RELEASE"
BUILD_DIR="./build"

# Set PETSc paths
PETSC_BUILD_DIR=$PETSC_DIR/$PETSC_ARCH

while getopts ":g" opt; do
  case ${opt} in
    g )
      BUILD_TYPE="DEBUG"
      BUILD_DIR="./build-debug"
      echo "Building in debug mode..."
      ;;
    * )
      echo "Usage: $0 [-g]"
      exit 1
      ;;
  esac
done

rm -rf $BUILD_DIR
cmake \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
  -DMFEM_USE_PETSC=YES \
  -DPETSC_DIR=$PETSC_DIR \
  -DPETSC_ARCH=$PETSC_ARCH \
  -DMFEM_USE_MPI=YES \
  -DHYPRE_DIR=$PETSC_BUILD_DIR \
  -DMFEM_USE_LAPACK=YES \
  -DBLAS_LIBRARIES=$PETSC_BUILD_DIR \
  -DLAPACK_LIBRARIES=$PETSC_BUILD_DIR \
  -DMFEM_USE_METIS_5=YES \
  -DMETIS_DIR=$PETSC_BUILD_DIR \
  -DParMETIS_DIR=$PETSC_BUILD_DIR \
  -DScaLAPACK_DIR=$PETSC_DIR/$PETSC_ARCH/externalpackages/git.scalapack/petsc-build \
  -DMFEM_USE_SUITESPARSE=YES \
  -DSuiteSparse_DIR=$PETSC_BUILD_DIR \
  -DMFEM_USE_MUMPS=YES \
  -DMUMPS_DIR=$PETSC_BUILD_DIR \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -S . -B $BUILD_DIR
cat > .clangd <<EOF
CompileFlags:
  CompilationDatabase: $BUILD_DIR
  Add:
    - -DMFEM_CONFIG_FILE="$BUILDDIR/config/_config.hpp"
EOF
cd $BUILD_DIR
make mfem -j8
