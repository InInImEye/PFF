//+
rad = DefineNumber[ 2.5, Name "Parameters/rad" ];
//+
buff = DefineNumber[ .5, Name "Parameters/buff" ];
//+
refB = DefineNumber[ 5., Name "Parameters/refB" ];
//+
refS = DefineNumber[ .2, Name "Parameters/refS" ];
//+
refs = DefineNumber[ .1, Name "Parameters/refs" ];
//+
Point(1) = {-9, -25, 0, refB};
//+
Point(2) = {9, -25, 0, refB};
//+
Point(3) = {9, 25, 0, refB};
//+
Point(4) = {-9, 25, 0, refB};
//+
Point(5) = {0, 0, 0, refs};
//+
Point(6) = {-9 + rad, 0, 0, refs};
//+
Point(7) = {9 - rad, 0, 0, refs};
//+
Point(8) = {9, rad, 0, refS};
//+
Point(9) = {9, -rad, 0, refS};
//+
Point(10) = {-9, -rad, 0, refS};
//+
Point(11) = {-9, rad, 0, refS};
//+
Point(12) = {0, rad, 0, refS};
//+
Point(13) = {0, -rad, 0, refS};
//+
Point(14) = {0, 25, 0, refB};
//+
Point(15) = {0, -25, 0, refB};
//+
Point(16) = {9, 0, 0, refB};
//+
Point(17) = {-9, 0, 0, refB};
//+
Line(1) = {1, 15};
//+
Line(2) = {15, 2};
//+
Line(3) = {2, 9};
//+
Circle(4) = {9, 16, 7};
//+
Circle(5) = {7, 16, 8};
//+
Line(6) = {8, 3};
//+
Line(7) = {3, 14};
//+
Line(8) = {14, 4};
//+
Line(9) = {4, 11};
//+
Circle(10) = {11, 17, 6};
//+
Circle(11) = {6, 17, 10};
//+
Line(12) = {10, 1};
//+
Line(13) = {15, 13};
//+
Line(14) = {13, 10};
//+
Line(15) = {9, 13};
//+
Line(16) = {13, 5};
//+
Line(17) = {5, 6};
//+
Line(18) = {7, 5};
//+
Line(19) = {5, 12};
//+
Line(20) = {12, 11};
//+
Line(21) = {8, 12};
//+
Line(22) = {12, 14};
//+
Curve Loop(1) = {1, 13, 14, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 15, -13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {14, -11, -17, -16};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {15, 16, -18, -4};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {18, 19, -21, -5};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {17, -10, -20, -19};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {20, -9, -8, -22};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {21, 22, -7, -6};
//+
Plane Surface(8) = {8};
//+
Physical Surface(2) = {8, 5, 4, 2, 1, 3, 6, 7};
//+
Physical Curve(0) = {3, 4, 5, 6, 9, 10, 11, 12};
//+
Physical Curve(11) = {7, 8};
//+
Physical Curve(13) = {1, 2};

// Meshing

Mesh.Algorithm = 3;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 1.;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
