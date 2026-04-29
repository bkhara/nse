GMSH="gmsh"

# geo=fpc2d_mesh_0.geo
geo=fpc_uns.geo

${GMSH} -format msh22 -2 "$geo" -o "./$(basename "${geo%.geo}").msh"
