// -----------------------------------------------------------------------------
// Lid-driven cavity mesh: unit square [0,1] x [0,1]
// Triangular transfinite mesh with clustering toward all four walls.
//
// Boundary physical tags match PCase_LDC_2D.h:
//   LEFT  = 1
//   RIGHT = 2
//   BTM   = 3
//   TOP   = 4
// Domain physical tag = 1
//
// Generate with, for example:
//   gmsh -2 ldc_unit_square_tri_wall_refined.geo -format msh2
//   meshio convert ldc_unit_square_tri_wall_refined.msh ldc_unit_square_tri_wall_refined.mesh
// -----------------------------------------------------------------------------

SetFactory("Built-in");

L = 1.0;

// Number of points per coordinate direction. 65 gives 64 intervals in x and y.
// Increase to 97 or 129 if you want a denser cavity mesh.
N = 65;

// Bump < 1 clusters points toward both ends of each curve.  Smaller values give
// stronger wall/corner clustering.  0.12 is a moderate boundary refinement.
bump = 0.12;

// Geometry
Point(1) = {0, 0, 0, 1.0};
Point(2) = {L, 0, 0, 1.0};
Point(3) = {L, L, 0, 1.0};
Point(4) = {0, L, 0, 1.0};

Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Transfinite grid with symmetric clustering at both ends of each coordinate
// direction.  Since the surface is not recombined, Gmsh splits the structured
// grid into triangles.
Transfinite Curve {1, 3} = N Using Bump bump;
Transfinite Curve {2, 4} = N Using Bump bump;
Transfinite Surface {1} = {1, 2, 3, 4};

// Triangles only.  Do not recombine into quads.
Mesh.RecombineAll = 0;
Recombine Surface {}; // intentionally empty; keeps triangles

// Physical tags expected by the LDC problem case.
Physical Curve("left",   1) = {4};
Physical Curve("right",  2) = {2};
Physical Curve("bottom", 3) = {1};
Physical Curve("top",    4) = {3};
Physical Surface("fluid", 1) = {1};

// Mesh output options.
Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 2.2;
