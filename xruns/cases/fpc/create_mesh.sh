GMSH="gmsh"

# geo=fpc_uns.geo
geo=fpc_str_0.geo
# geo=fpc_str_1.geo

${GMSH} -format msh22 -2 "$geo" -o "./$(basename "${geo%.geo}").msh"
