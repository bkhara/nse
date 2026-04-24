GMSH="gmsh"

geo=fpc2d_mesh_0.geo

${GMSH} -format msh22 -2 "$geo" -o "./$(basename "${geo%.geo}").msh"
