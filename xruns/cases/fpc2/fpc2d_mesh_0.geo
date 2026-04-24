// PARAMTER INITIALIZATIONS
_centerX = 0;
_centerY = 0;
_radius = 0.5;

_leftWallFromCenter = 2.5;
_rightWallFromCenter = 30;

_partLineFromCenter = 2.5;

_bottomWallFromCenter = 2.5;
_topWallFromCenter = 2.5;

_N_circle = 21;
_R_circle = 1.00;

_N_partition_lines = 31;
_R_partition_lines = 0.97;

_N_wake = 161;
_R_wake = 1.00;

_zHeight = 0.1;
_zNumLayers = 1;


// GEOMETRIC PARAMETERS

centerX 	= DefineNumber[ _centerX, 	Name "Parameters/Geometry/circle/centerX" ];
centerY	    = DefineNumber[ _centerY, 	Name "Parameters/Geometry/circle/centerY" ];
radius	    = DefineNumber[ _radius, 	Name "Parameters/Geometry/circle/radius" ];
leftWallFromCenter	    = DefineNumber[ _leftWallFromCenter, 	Name "Parameters/Geometry/Wall/leftWallFromCenter" ];
rightWallFromCenter	    = DefineNumber[ _rightWallFromCenter, 	Name "Parameters/Geometry/Wall/rightWallFromCenter" ];
bottomWallFromCenter	= DefineNumber[ _bottomWallFromCenter, 	Name "Parameters/Geometry/Wall/bottomWallFromCenter" ];
topWallFromCenter	    = DefineNumber[ _topWallFromCenter, 	Name "Parameters/Geometry/Wall/topWallFromCenter" ];
partLineFromCenter	    = DefineNumber[ _partLineFromCenter, 	Name "Parameters/Geometry/PartitionLine/partLineFromCenter" ];

// MESH PARAMETERS

N_circle = DefineNumber[ _N_circle, 	Name "Parameters/Mesh/Circle/N_circle" ];
R_circle = DefineNumber[ _R_circle, 	Name "Parameters/Mesh/Circle/R_circle" ];

N_partition_lines 	= DefineNumber[ _N_partition_lines, 	Name "Parameters/Mesh/PartitionLines/N_partition_lines" ];
R_partition_lines 	= DefineNumber[ _R_partition_lines, 	Name "Parameters/Mesh/PartitionLines/R_partition_lines" ];

N_wake = DefineNumber[ _N_wake, 	Name "Parameters/Mesh/Wake/N_wake" ];
R_wake = DefineNumber[ _R_wake, 	Name "Parameters/Mesh/Wake/R_wake" ];

zHeight    = DefineNumber[ _zHeight, 	Name "Parameters/Mesh/Extrusion/zHeight" ];
zNumLayers = DefineNumber[ _zNumLayers, Name "Parameters/Mesh/Extrusion/zNumLayers" ];

theta = Atan(leftWallFromCenter / bottomWallFromCenter);
del_x = radius * Sin(theta);
del_y = radius * Cos(theta);

//rectangle
Point(1) = {centerX - leftWallFromCenter,  centerY - bottomWallFromCenter, 0, 1.0};
Point(2) = {centerX + partLineFromCenter,  centerY - bottomWallFromCenter, 0, 1.0};
Point(3) = {centerX + rightWallFromCenter, centerY - bottomWallFromCenter, 0, 1.0};
Point(4) = {centerX - leftWallFromCenter,  centerY + topWallFromCenter, 0, 1.0};
Point(5) = {centerX + partLineFromCenter,  centerY + topWallFromCenter, 0, 1.0};
Point(6) = {centerX + rightWallFromCenter, centerY + topWallFromCenter, 0, 1.0};

// circle
Point(7)  = {centerX - del_x, centerY - del_y, 0, 1.0};
Point(8)  = {centerX + del_x, centerY - del_y, 0, 1.0};
Point(9)  = {centerX, centerY, 0, 1.0};
Point(10) = {centerX - del_x, centerY + del_y, 0, 1.0};
Point(11) = {centerX + del_x, centerY + del_y, 0, 1.0};

// rectangle lines
Line(1) = {1, 2}; Transfinite Line {1} = N_circle Using Bump R_circle; // bottom left line - horizontal
Line(3) = {4, 5}; Transfinite Line {3} = N_circle Using Bump R_circle; // top left line - horizontal

Line(2) = {2, 3}; Transfinite Line {2} = N_wake Using Bump R_wake; // bottom right line - horizontal
Line(4) = {5, 6}; Transfinite Line {4} = N_wake Using Bump R_wake; // top right line - horizontal

Line(5) = {1, 4}; Transfinite Line {5} = N_circle Using Bump R_circle; // left line - vertical
Line(6) = {2, 5}; Transfinite Line {6} = N_circle Using Bump R_circle; // middle partition line - vertical
Line(7) = {3, 6}; Transfinite Line {7} = N_circle Using Bump R_circle; // right line - vertical

Circle(8) = {7, 9, 8};    Transfinite Line {8}  = N_circle Using Bump R_circle; // bottom arc
Circle(9) = {8, 9, 11};   Transfinite Line {9}  = N_circle Using Bump R_circle; // right arc
Circle(10) = {11, 9, 10}; Transfinite Line {10} = N_circle Using Bump R_circle; // top arc
Circle(11) = {10, 9, 7};  Transfinite Line {11} = N_circle Using Bump R_circle; // left arc

Line(12) = {1, 7};  Transfinite Line {12} = N_partition_lines Using Progression R_partition_lines; // partition around circle - bottom left diagonal
Line(13) = {2, 8};  Transfinite Line {13} = N_partition_lines Using Progression R_partition_lines; // partition around circle - bottom right diagonal
Line(14) = {5, 11}; Transfinite Line {14} = N_partition_lines Using Progression R_partition_lines; // partition around circle - top right diagonal
Line(15) = {4, 10}; Transfinite Line {15} = N_partition_lines Using Progression R_partition_lines; // partition around circle - top left diagonal

// surfaces
Line Loop(16) = {12, 8, -13, -1};
Plane Surface(17) = {16};
Line Loop(18) = {13, 9, -14, -6};
Plane Surface(19) = {18};
Line Loop(20) = {14, 10, -15, 3};
Plane Surface(21) = {20};
Line Loop(22) = {15, 11, -12, 5};
Plane Surface(23) = {22};
Line Loop(24) = {2, 7, -4, -6};
Plane Surface(25) = {24};

Transfinite Surface {17};
Transfinite Surface {19};
Transfinite Surface {21};
Transfinite Surface {23};
Transfinite Surface {25};
Recombine Surface {17};
Recombine Surface {19};
Recombine Surface {21};
Recombine Surface {23};
Recombine Surface {25};


Physical Line(1) = {5};
Physical Line(2) = {7};
Physical Line(3) = {1, 2};
Physical Line(4) = {3, 4};
Physical Line(7) = {8, 9, 10, 11};
Physical Surface(1) = {17, 19, 21, 23, 25};
