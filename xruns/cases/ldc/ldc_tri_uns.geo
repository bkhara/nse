// -----------------------------------------------------------------------------
// Lid-driven cavity mesh: unit square [0,1] x [0,1]
// Unstructured triangular mesh with refinement near all boundaries.
//
// Boundary physical tags match PCase_LDC_2D.h:
//   LEFT   = 1
//   RIGHT  = 2
//   BTM    = 3
//   TOP    = 4
// Domain physical tag:
//   FLUID  = 1
// -----------------------------------------------------------------------------

SetFactory("OpenCASCADE");

// ----------------------------
// Geometry
// ----------------------------
L = 1.0;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {L,   0.0, 0.0};
Point(3) = {L,   L,   0.0};
Point(4) = {0.0, L,   0.0};

Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// ----------------------------
// Physical tags
// ----------------------------
Physical Curve("left",   1) = {4};
Physical Curve("right",  2) = {2};
Physical Curve("bottom", 3) = {1};
Physical Curve("top",    4) = {3};
Physical Surface("fluid", 1) = {1};

// ----------------------------
// Mesh-size field
// ----------------------------
// h_wall: small elements near all cavity walls
// h_bulk: larger elements near the center
// d_min : distance from wall where refinement is strongest
// d_max : distance from wall where mesh reaches h_bulk
//
// Tune these four values as needed.
h_wall = 0.003;
h_bulk = 0.01;
d_min  = 0.015;
d_max  = 0.25;

// Compute distance to all four wall curves.
Field[1] = Distance;
Field[1].CurvesList = {1, 2, 3, 4};
Field[1].Sampling = 200;

// Convert distance-to-wall into target element size.
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_wall;
Field[2].SizeMax = h_bulk;
Field[2].DistMin = d_min;
Field[2].DistMax = d_max;

Background Field = 2;

// Make the background field control the sizes instead of point/curvature defaults.
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Unstructured triangles.
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombineAll = 0;
Mesh.ElementOrder = 1;

// Save as Gmsh v2 ASCII if using `gmsh -2 -format msh2`.
Mesh.MshFileVersion = 2.2;
