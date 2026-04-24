./configure \
  --with-debugging=0 \
  --download-hypre \
  --download-mumps \
  --download-suitesparse \
  --download-openblas \
  --download-scalapack \
  --download-fblaslapack \
  --download-metis \
  --download-parmetis \
  FOPTFLAGS="-O3" \
  COPTFLAGS="-O3 -march=native -mtune=native" \
  CXXOPTFLAGS="-O3 -march=native -mtune=native"