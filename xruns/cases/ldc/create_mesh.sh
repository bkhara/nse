GMSH="gmsh"

geo=ldc_quad_bump.geo
#geo=ldc_tri_str.geo
#geo=ldc_tri_uns.geo

${GMSH} -format msh22 -2 "$geo" -o "./$(basename "${geo%.geo}").msh"
