//+
refs = DefineNumber[ 0.5, Name "Parameters/refs" ];
//+
refS = DefineNumber[ 0.7, Name "Parameters/refS" ];
//+
refB = DefineNumber[ 10., Name "Parameters/refB" ];
//+
Point(1) = {-30, -60, 0, refB};
//+
Point(2) = {30, -60, 0, refB};
//+
Point(3) = {30, 60, 0, refB};
//+
Point(4) = {-30, 60, 0, refB};
//+
Point(5) = {0, 60, 0, refB};
//+
Point(6) = {0, -60, 0, refB};
//+
Point(7) = {0, 0, 0, refs};
//+
Point(8) = {0, 10, 0, refS};
//+
Point(9) = {0, -10, 0, refS};
//+
Point(10) = {-30, -40, 0, refB};
//+
Point(11) = {-30, -30, 0, refB};
//+
Point(12) = {-20, -30, 0, refB};
//+
Point(13) = {20, -30, 0, refB};
//+
Point(14) = {30, -30, 0, refB};
//+
Point(15) = {30, -40, 0, refB};
//+
Point(16) = {30, 40, 0, refB};
//+
Point(17) = {30, 30, 0, refB};
//+
Point(18) = {20, 30, 0, refB};
//+
Point(19) = {-20, 30, 0, refB};
//+
Point(20) = {-30, 30, 0, refB};
//+
Point(21) = {-30, 40, 0, refB};
//+
Point(22) = {-20, 0, 0, refs};
//+
Point(23) = {20, 0, 0, refs};
//+
Point(24) = {20, 10, 0, refS};
//+
Point(25) = {-20, 10, 0, refS};
//+
Point(26) = {-20, -10, 0, refS};
//+
Point(27) = {20, -10, 0, refS};
//+
Point(28) = {0, -30, 0, refB};
//+
Point(29) = {0, 30, 0, refB};
//+
Line(1) = {1, 6};
//+
Line(2) = {6, 2};
//+
Line(3) = {2, 15};
//+
Circle(4) = {15, 14, 13};
//+
Line(5) = {13, 27};
//+
Line(6) = {27, 23};
//+
Line(7) = {23, 24};
//+
Line(8) = {24, 18};
//+
Circle(9) = {18, 17, 16};
//+
Line(10) = {16, 3};
//+
Line(11) = {3, 5};
//+
Line(12) = {5, 4};
//+
Line(13) = {4, 21};
//+
Circle(14) = {21, 20, 19};
//+
Line(15) = {19, 25};
//+
Line(16) = {25, 22};
//+
Line(17) = {22, 26};
//+
Line(18) = {26, 12};
//+
Circle(19) = {12, 11, 10};
//+
Line(20) = {10, 1};
//+
Line(21) = {6, 28};
//+
Line(22) = {28, 12};
//+
Line(23) = {13, 28};
//+
Line(24) = {28, 9};
//+
Line(25) = {9, 26};
//+
Line(26) = {27, 9};
//+
Line(27) = {9, 7};
//+
Line(28) = {7, 22};
//+
Line(29) = {23, 7};
//+
Line(30) = {7, 8};
//+
Line(31) = {8, 25};
//+
Line(32) = {24, 8};
//+
Line(33) = {8, 29};
//+
Line(34) = {29, 19};
//+
Line(35) = {18, 29};
//+
Line(36) = {29, 5};
//+
Curve Loop(1) = {1, 21, 22, 19, 20};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 4, 23, -21};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {22, -18, -25, -24};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {23, 24, -26, -5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {25, -17, -28, -27};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {26, 27, -29, -6};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {28, -16, -31, -30};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {29, 30, -32, -7};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {31, -15, -34, -33};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {32, 33, -35, -8};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {34, -14, -13, -12, -36};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {35, 36, -11, -10, -9};
//+
Plane Surface(12) = {12};
//+
Physical Surface(2) = {12, 10, 8, 6, 4, 2, 1, 3, 5, 7, 9, 11};
//+
Physical Curve(0) = {3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 19, 20};
//+
Physical Curve(11) = {11, 12};
//+
Physical Curve(13) = {1, 2};

// Meshing

Mesh.Algorithm = 3;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 1.;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
