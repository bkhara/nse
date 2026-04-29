// -----------------------------------------------------------------------------
//
//  Gmsh GEO tutorial 18
//
//  Periodic meshes
//
// -----------------------------------------------------------------------------

// Periodic meshing constraints can be imposed on surfaces and curves.

// Let's use the OpenCASCADE geometry kernel to build two geometries.

SetFactory("OpenCASCADE");

// -----------------------------------------------------------------------------
// Geometry parameters
// -----------------------------------------------------------------------------

C = 1;

Origin_X = 0;
Origin_Y = 0;
Origin_Z = 0;

Lx = 20;
Ly = 5;
Lz = 0;

// -----------------------------------------------------------------------------
// Symmetry-breaking parameter
//
// If circle_y_shift = 0, the geometry is symmetric about y = 0.
// Positive values move the cylinder upward.
// Negative values move the cylinder downward.
// -----------------------------------------------------------------------------

circle_y_shift = -0.1;

// Channel placement.
// With CY = Ly/2 and ylo = Origin_Y - CY,
// the channel extends from y = -Ly/2 to y = +Ly/2.
CX = 3;       // center of circle measured from x-low of channel
CY = Ly/2;    // half channel height
CZ = Lz/2;

xlo = Origin_X - CX;
ylo = Origin_Y - CY;
zlo = Origin_Z - CZ;

// Circle center.
// If circle_y_shift = 0, this is exactly (Origin_X, Origin_Y, Origin_Z).
Circle_CX = Origin_X;
Circle_CY = Origin_Y + circle_y_shift;
Circle_CZ = Origin_Z;

// Mesh sizes
lc = 0.5;
lc_af = 0.5;

// -----------------------------------------------------------------------------
// Geometry
// -----------------------------------------------------------------------------

Rectangle(1) = {xlo, ylo, 0, Lx, Ly, 0};

// Obstacle/cylinder
Sradius = 0.5;

Circle(7) = {Circle_CX, Circle_CY, Circle_CZ, Sradius, 0, 2*Pi};

Curve Loop(2) = {7};

Plane Surface(2) = {2};

// Remove cylinder from channel region
s = 9;
BooleanDifference(s) = { Surface{1}; Delete; }{ Surface{2}; Delete; };

// -----------------------------------------------------------------------------
// Physical groups
//
// Configure which curves/surfaces get boundary/domain attributes.
// These curve IDs are assumed from the current BooleanDifference result.
// If you change the geometry substantially, verify with Gmsh GUI.
// -----------------------------------------------------------------------------

Physical Line(1) = {2}; // -x
Physical Line(2) = {3}; // +x
Physical Line(3) = {1}; // -y
Physical Line(4) = {4}; // +y
Physical Curve(7) = {5}; // cylinder/stator

Physical Surface(1) = {s};

// -----------------------------------------------------------------------------
// Mesh configuration
// -----------------------------------------------------------------------------

// Refinement around the cylinder boundary
F1 = 1;
Field[F1] = Distance;
Field[F1].CurvesList = {5};
Field[F1].Sampling = 100;

F2 = 2;
Field[F2] = Threshold;
Field[F2].InField = F1;
Field[F2].SizeMin = lc / 10;
Field[F2].SizeMax = lc / 2;
Field[F2].DistMin = 1.0 * C;
Field[F2].DistMax = 2.0 * C;

// Refinement along a line through the cylinder center.
// If circle_y_shift = 0, this line is y = 0.
// If circle_y_shift != 0, the refinement band follows the shifted cylinder.
p = newp;
nDistPts = 50;

For i In {0:nDistPts}
     x = xlo + i * Lx / nDistPts;
     Point(p + i) = {x, circle_y_shift, 0, lc_af};
EndFor

// Rotate {{0, 0, -1}, {0.5, 0, 0}, AoA} { Point{p:p+nDistPts}; }

F6 = 6;
Field[F6] = Distance;
Field[F6].PointsList = {p:p+nDistPts};

F7 = 7;
Field[F7] = Threshold;
Field[F7].InField = F6;
Field[F7].SizeMin = lc / 10;
Field[F7].SizeMax = lc / 4;
Field[F7].DistMin = 1.0 * C;
Field[F7].DistMax = 2.0 * C;

// -----------------------------------------------------------------------------
// Optional box refinement fields, currently disabled
// -----------------------------------------------------------------------------

// Field[1] = Box;
// Field[1].VIn = lc/5;
// Field[1].VOut = lc/2;
// Field[1].XMin = origin_x + 0;
// Field[1].XMax = origin_x + Lx;
// Field[1].YMin = - Cr/3;
// Field[1].YMax = + Cr/3;
// Field[1].ZMin = 0;
// Field[1].ZMax = 0;

// Field[2] = Box;
// Field[2].VIn  = lc/10;
// Field[2].VOut = lc;
// Field[2].XMin = airfoil_offset_x - Cr/10;
// Field[2].XMax = airfoil_offset_x + Cr + Cr/10;
// Field[2].YMin = airfoil_offset_y - Cr/10;
// Field[2].YMax = airfoil_offset_y + Cr/6;
// Field[2].ZMin = 0;
// Field[2].ZMax = 0;

// -----------------------------------------------------------------------------
// Background mesh field
// -----------------------------------------------------------------------------

Field[10] = Min;
Field[10].FieldsList = {F2, F7};

Background Field = 10;
