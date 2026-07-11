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

# Available methods

A "method" is a combination of four independent choices:

* **coupling** — `coupled` (monolithic velocity-pressure) or `uncoupled` (projection / Chorin-Temam).
* **marching** — `bdf1` or `bdf2` (default `bdf2`). Crank-Nicolson (`cn`) is retired: the CN
  integrators remain in the source as dead code but are no longer selectable.
* **convection form** — `conv` (advective), `skew` (skew-symmetric), or `div` (divergence / conservative).
* **stabilization** — `none`, `vms` (residual-based VMS), or `sups` (SUPG + PSPG).

Only a subset of the full cross-product is wired up. The table below is the single source of
truth for what actually runs; anything not listed is rejected at startup.

| Coupling   | Marching | Convection | Stabilization | Status      | Notes                                             |
|------------|----------|------------|---------------|-------------|---------------------------------------------------|
| coupled    | bdf2     | div        | none          | available   | div-form Galerkin, no stabilization               |
| coupled    | bdf2     | div        | vms           | available   | residual-based VMS                                |
| coupled    | bdf2     | div        | sups          | available   | SUPG + PSPG add-ons (see `use_supg` / `use_pspg`) |
| uncoupled  | bdf2     | conv       | vms           | available   | incremental-pressure projection                   |
| uncoupled  | bdf2     | div        | vms           | available   | incremental-pressure projection (+ outflow flux)  |
| uncoupled  | bdf1     | conv       | vms           | planned     | Chorin first-order projection                     |
| uncoupled  | bdf1     | div        | vms           | planned     | Chorin first-order projection                     |

Notes:

* **Coupled solver is fixed to BDF2 for now.** `NSEBlockIntegBDF2` hardcodes the BDF2 coefficients, so
  there is no coupled BDF1 path yet.
* **Projection temporal order follows `marching`.** `bdf1` maps to the Chorin first-order scheme and
  `bdf2` to the incremental-pressure BDF2 scheme; `projection_config.scheme` is derived from
  `marching` rather than set independently.
* **Projection is always VMS-stabilized.** The projection integrators carry the residual-based
  `tau_M` terms intrinsically, so `stab = none` / `sups` are not valid for the uncoupled path.
* **`use_supg` / `use_pspg`** are sub-toggles within `stab = sups`. Turning both off yields an
  unstabilized run (equivalent to `stab = none`).
* The stabilization parameter constants (`Ci`, etc.) are common to every stabilized method and are
  not part of the method selection.