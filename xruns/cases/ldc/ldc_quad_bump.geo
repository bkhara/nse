// -----------------------------------------------------------------------------
// Lid-driven cavity mesh: unit square [0,1] x [0,1]
// Structured quadrilateral mesh with bump-type boundary refinement.
//
// Physical tags:
//   left   = 1
//   right  = 2
//   bottom = 3
//   top    = 4
//   fluid  = 1
//
// Generate with:
//   gmsh -2 ldc_unit_square_quad_bump.geo -format msh2
// -----------------------------------------------------------------------------

SetFactory("Built-in");

// -----------------------------------------------------------------------------
// User controls
// -----------------------------------------------------------------------------
Nx = 128;        // number of quad elements in x
Ny = 128;        // number of quad elements in y
bump = 0.08;    // smaller => stronger clustering near both ends of each curve

Lx = 1.0;
Ly = 1.0;

// -----------------------------------------------------------------------------
// Geometry
// -----------------------------------------------------------------------------
Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {Lx,  0.0, 0.0, 1.0};
Point(3) = {Lx,  Ly,  0.0, 1.0};
Point(4) = {0.0, Ly,  0.0, 1.0};

Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// -----------------------------------------------------------------------------
// Structured quad mesh with wall refinement
//
// Using Bump clusters points toward both ends of a curve. Since the surface is
// transfinite, clustering the horizontal curves gives refinement near x=0 and
// x=1, while clustering the vertical curves gives refinement near y=0 and y=1.
// -----------------------------------------------------------------------------
Transfinite Curve {1, 3} = Nx + 1 Using Bump bump;
Transfinite Curve {2, 4} = Ny + 1 Using Bump bump;
Transfinite Surface {1} = {1, 2, 3, 4};
Recombine Surface {1};

// Optional smoothing improves quad quality slightly without changing topology.
Mesh.Smoothing = 5;

// -----------------------------------------------------------------------------
// Physical tags for MFEM/problem case compatibility
// -----------------------------------------------------------------------------
Physical Curve("left",   1) = {4};
Physical Curve("right",  2) = {2};
Physical Curve("bottom", 3) = {1};
Physical Curve("top",    4) = {3};
Physical Surface("fluid", 1) = {1};
