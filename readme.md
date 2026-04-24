# Building dependencies
## PETSc setup
* Installation
```
./configure --download-metis --download-parmetis --download-scalapack --download-fblaslapack --download-mumps --download-hypre -with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native'
# optional: --download-mpich
make
make install
```
* Export the following variables
```
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=<petsc-arch>
```
* Or use the script at [install_scripts/petsc.sh](install_scripts/petsc.sh)

## Libconfig setup
* Installation
```
wget https://hyperrealm.github.io/libconfig/dist/libconfig-1.7.3.tar.gz
tar xvf libconfig-1.7.3.tar.gz
rm libconfig-1.7.3.tar.gz
cd libconfig-1.7.3
./configure --prefix=`pwd`/install
make
make install
```
* Export the following variables
```
export LIBCONFIG_DIR=<path-to-libconfig>/install
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LIBCONFIG_DIR}/lib
```
# Building MFEM
```
~> ls
glvis-4.2/  hypre-2.26.0.tar.gz   metis-4.0.3.tar.gz   mfem-4.5/
```
### HYPRE
* HYPRE installation is usual:
```
~> tar -zxvf hypre-2.26.0.tar.gz
~> cd hypre-2.26.0/src/
~/hypre-2.26.0/src> ./configure --disable-fortran
~/hypre-2.26.0/src> make -j
~/hypre-2.26.0/src> cd ../..
~> ln -s hypre-2.26.0 hypre
```

### METIS-5
* METIS-5 installation is usual:
```
~> tar zvxf metis-5.1.0.tar.gz
~> cd metis-5.1.0
~/metis-5.1.0> make BUILDDIR=lib config
~/metis-5.1.0> make BUILDDIR=lib
~/metis-5.1.0> cp lib/libmetis/libmetis.a lib
```

### MFEM installation
* Method I 
```
~> cd mfem-4.5
~/mfem-4.5> make parallel -j MFEM_USE_METIS_5=YES METIS_DIR=@MFEM_DIR@/../metis-5.1.0
```
* Method II (using CMake)
```
cd mfem-4.8
mkdir build; cd build
cmake .. -DMFEM_USE_MPI=ON -DMFEM_USE_METIS_5=YES -DMETIS_DIR=/Users/biswajitkhara/packages/mfem/metis-5.1.0
cmake .. -DMFEM_USE_MPI=ON -DMFEM_USE_METIS_5=YES -DMFEM_USE_PETSC=YES -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH -DMETIS_DIR=$METIS_DIR
cmake .. -DMFEM_USE_MPI=ON -DMFEM_USE_METIS_5=YES -DMFEM_USE_PETSC=YES -DPETSC_DIR=/Users/biswajitkhara/packages/petsc-3.23.4 -DPETSC_ARCH=arch-darwin-c-opt -DMETIS_DIR=/Users/biswajitkhara/packages/mfem/metis-5.1.0
make -j8
```
* Or use the script at [install_scripts/mfem.sh](install_scripts/mfem.sh)

# Building this code (fracture)
```
mkdir build; cd build
cmake ..
make
```