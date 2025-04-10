//+
buff = DefineNumber[ 0.5, Name "Parameters/buff" ];
//+
refB = DefineNumber[ 4., Name "Parameters/refB" ];
//+
refS = DefineNumber[ .2, Name "Parameters/refS" ];
//+
refs = DefineNumber[ .05, Name "Parameters/refs" ];
//+
rad = DefineNumber[ 2.5, Name "Parameters/rad" ];
//+
in = DefineNumber[ 1., Name "Parameters/in" ];
//+
thet = DefineNumber[ Acos(in / rad), Name "Parameters/thet" ];
//+
len = DefineNumber[ 5., Name "Parameters/len" ];

//+
Point(1) = {0., 0., 0., refs};
//+
coord = len - (rad * Cos(Pi / 4) + in);
//+
Point(2) = {coord, coord, 0, refs};
//+
Point(3) = {-coord, -coord, 0, refs};
//+
curend = len -(in + rad * Sin(thet));
//+
Point(4) = {curend, len, 0, refS};
//+
Point(5) = {len, curend, 0, refS};
//+
Point(6) = {-len, -curend, 0, refS};
//+
Point(7) = {-curend, -len, 0, refS};
//+
Point(8) = {len - in, len -in, 0, refs};
//+
Point(9) = {in - len, in - len, 0, refs};
//+
Point(10) = {curend-buff, len, 0, refS};
//+
Point(11) = {len, curend-buff, 0, refS};
//+
Point(12) = {-len, -curend+buff, 0, refS};
//+
Point(13) = {-curend+buff, -len, 0, refS};
//+
Point(14) = {-len, len, 0, refB};
//+
Point(15) = {len, -len, 0, refB};
//+
Line(1) = {7, 13};
//+
Line(2) = {13, 15};
//+
Line(3) = {15, 11};
//+
Line(4) = {11, 5};
//+
Circle(5) = {5, 8, 2};
//+
Circle(6) = {2, 8, 4};
//+
Line(7) = {4, 10};
//+
Line(8) = {10, 14};
//+
Line(9) = {14, 12};
//+
Line(10) = {12, 6};
//+
Circle(11) = {6, 9, 3};
//+
Circle(12) = {3, 9, 7};
//+
Line(13) = {3, 1};
//+
Line(14) = {1, 2};
//+
Line(15) = {12, 10};
//+
Line(16) = {13, 11};
//+
Line(17) = {1, 10};
//+
Line(18) = {1, 12};
//+
Line(19) = {1, 13};
//+
Line(20) = {1, 11};
//+
Curve Loop(1) = {1, -19, -13, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {16, -20, 19};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, 3, -16};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {20, 4, 5, -14};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {11, 13, 18, 10};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {18, 15, -17};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {14, 6, 7, -17};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {15, 8, 9};
//+
Plane Surface(8) = {8};
//+
Physical Surface(2) = {3, 2, 1, 5, 6, 8, 7, 4};
//+
Physical Curve(0) = {11, 12, 5, 6};
//+
Physical Curve(13) = {1, 2, 3, 4};
//+
Physical Curve(11) = {7, 8};
//+
Physical Curve(12) = {9, 10};


// Meshing

Mesh.Algorithm = 3;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 1.;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
