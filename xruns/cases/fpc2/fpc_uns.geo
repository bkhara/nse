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

// The first geometry is very simple: a unit cube with a non-uniform mesh size
// constraint (set on purpose to be able to verify visually that the periodicity
// constraint works!):
C = 1;
Origin_X = 0;
Origin_Y = 0;
Origin_Z = 0;

Lx = 20;
Ly = 5;
Lz = 0;

CX = 3; // center of sphere from Xlow of channel
CY = Ly/2; // center of sphere from Ylow of channel
CZ = Lz/2; // center of sphere from Zlow of channel
xlo = Origin_X - CX;
ylo = Origin_Y - CY;
zlo = Origin_Z - CZ;

// Printf("CX = %f", CX);
// Printf("CY = %f", CY);
// Printf("CZ = %f", CZ);
// Printf("xlo = %f", xlo);
// Printf("ylo = %f", ylo);
// Printf("zlo = %f", zlo);

lc = 0.5;
lc_af = 0.5;

Rectangle(1) = {xlo, ylo, 0, Lx, Ly, 0};

// We start with a cube and some spheres:
Sradius = 0.5;
//+
Circle(7) = {Origin_X, Origin_Y, Origin_Z, Sradius, 0, 2*Pi};

//+
Curve Loop(2) = {7};
//+
Plane Surface(2) = {2};

s = 9;
BooleanDifference(s) = { Surface{1}; Delete; }{ Surface{2}; Delete; }; //remove rotor from the blockmesh region

//Configure for which surfaces we should generate face elements (and in which order),
//without it TalyFEM gives error on mesh loading
Physical Line(1) = {2}; //-x
Physical Line(2) = {3}; //+x
Physical Line(3) = {1};  //-y
Physical Line(4) = {4};  //+y
Physical Curve(7) = {5};  //stators

Physical Surface(1) = {s};

// mesh configuration
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

p = newp;
nDistPts = 50;
For i In {0: nDistPts}
     x = xlo + i * Lx / nDistPts;
     Point(p + i) = {x, 0, 0, lc_af};
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
// Field[2].YMin = airfoil_offset_y-Cr/10;
// Field[2].YMax = airfoil_offset_y+Cr/6;
// Field[2].ZMin = 0;
// Field[2].ZMax = 0;

Field[10] = Min;
// Field[10].FieldsList = {1};
Field[10].FieldsList = {F2,F7};
Background Field = 10;