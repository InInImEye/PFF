//+
buff = DefineNumber[ .5, Name "Parameters/buff" ];
//+
refB = DefineNumber[ 8., Name "Parameters/refB" ];
//+
refS = DefineNumber[ .4, Name "Parameters/refS" ];
//+
refs = DefineNumber[ .1, Name "Parameters/refs" ];
//+
Point(1) = {-4, -50, 0, refB};
//+
Point(2) = {-4, 50, 0, refB};
//+
Point(3) = {4, 50, 0, refB};
//+
Point(4) = {4, -50, 0, refB};
//+
Point(5) = {1, 1, 0, refs};
//+
Point(6) = {1, -1, 0, refs};
//+
Point(7) = {-1, -1, 0, refs};
//+
Point(8) = {-1, 1, 0, refs};
//+
Point(9) = {1.5, 1.5, 0, refs};
//+
Point(10) = {1.5, -1.5, 0, refs};
//+
Point(11) = {-1.5, -1.5, 0, refs};
//+
Point(12) = {-1.5, 1.5, 0, refs};
//+
Point(13) = {1.5, 1., 0, refs};
//+
Point(14) = {1.5, -1., 0, refs};
//+
Point(15) = {-1.5, -1., 0, refs};
//+
Point(16) = {-1.5, 1., 0, refs};
//+
Point(17) = {0, 0, 0, refs};
//+
Point(18) = {1, 0, 0, refs};
//+
Point(19) = {-1, 0, 0, refs};
//+
Point(20) = {0, 1, 0, refs};
//+
Point(21) = {0, -1, 0, refs};
//+
Point(22) = {0, 2.5, 0, refs};
//+
Point(23) = {0, -2.5, 0, refs};
//+
Point(24) = {1.5, -2.5, 0, refs};
//+
Point(25) = {-1.5, -2.5, 0, refs};
//+
Point(26) = {-1.5, 2.5, 0, refs};
//+
Point(27) = {1.5, 2.5, 0, refs};
//+
Point(28) = {4, 1.5, 0, refS};
//+
Point(29) = {-4, 1.5, 0, refS};
//+
Point(30) = {-4, -1.5, 0, refS};
//+
Point(31) = {4, -1.5, 0, refS};
//+
Point(32) = {4, 2.5, 0, refS};
//+
Point(33) = {-4, 2.5, 0, refS};
//+
Point(34) = {-4, -2.5, 0, refS};
//+
Point(35) = {4, -2.5, 0, refS};
//+
Point(36) = {4, 3.5, 0, buff};
//+
Point(37) = {4, -3.5, 0, buff};
//+
Point(38) = {-4, -3.5, 0, buff};
//+
Point(39) = {-4, 3.5, 0, buff};
//+
Point(40) = {0, 3.5, 0, buff};
//+
Point(41) = {0, -3.5, 0, buff};
//+
Point(42) = {0, 50, 0, refB};
//+
Point(43) = {0, -50, 0, refB};
//+
Line(1) = {1, 43};
//+
Line(2) = {43, 4};
//+
Line(3) = {4, 37};
//+
Line(4) = {37, 35};
//+
Line(5) = {35, 31};
//+
Line(6) = {31, 10};
//+
Circle(7) = {10, 14, 6};
//+
Line(8) = {6, 18};
//+
Line(9) = {18, 5};
//+
Circle(10) = {5, 13, 9};
//+
Line(11) = {9, 28};
//+
Line(12) = {28, 32};
//+
Line(13) = {32, 36};
//+
Line(14) = {36, 3};
//+
Line(15) = {3, 42};
//+
Line(16) = {42, 2};
//+
Line(17) = {2, 39};
//+
Line(18) = {39, 33};
//+
Line(19) = {33, 29};
//+
Line(20) = {29, 12};
//+
Circle(21) = {12, 16, 8};
//+
Line(22) = {8, 19};
//+
Line(23) = {19, 7};
//+
Circle(24) = {7, 15, 11};
//+
Line(25) = {11, 30};
//+
Line(26) = {30, 34};
//+
Line(27) = {34, 38};
//+
Line(28) = {38, 1};
//+
Line(29) = {43, 41};
//+
Line(30) = {41, 38};
//+
Line(31) = {37, 41};
//+
Line(32) = {41, 23};
//+
Line(33) = {23, 25};
//+
Line(34) = {25, 34};
//+
Line(35) = {35, 24};
//+
Line(36) = {24, 23};
//+
Line(37) = {25, 11};
//+
Line(38) = {23, 21};
//+
Line(39) = {21, 7};
//+
Line(40) = {24, 10};
//+
Line(41) = {6, 21};
//+
Line(42) = {21, 17};
//+
Line(43) = {17, 19};
//+
Line(44) = {18, 17};
//+
Line(45) = {17, 20};
//+
Line(46) = {20, 8};
//+
Line(47) = {5, 20};
//+
Line(48) = {20, 22};
//+
Line(49) = {22, 26};
//+
Line(50) = {26, 12};
//+
Line(51) = {26, 33};
//+
Line(52) = {9, 27};
//+
Line(53) = {27, 22};
//+
Line(54) = {32, 27};
//+
Line(55) = {22, 40};
//+
Line(56) = {40, 39};
//+
Line(57) = {36, 40};
//+
Line(58) = {40, 42};
//+
Curve Loop(1) = {1, 29, 30, 28};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 31, -29};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {30, -27, -34, -33, -32};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {31, 32, -36, -35, -4};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {34, -26, -25, -37};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {33, 37, -24, -39, -38};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {36, 38, -41, -7, -40};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {35, 40, -6, -5};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {39, -23, -43, -42};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {41, 42, -44, -8};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {43, -22, -46, -45};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {44, 45, -47, -9};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {46, -21, -50, -49, -48};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {47, 48, -53, -52, -10};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {20, -50, 51, 19};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {11, 12, 54, -52};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {51, -18, -56, -55, 49};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {53, 55, -57, -13, 54};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {56, -17, -16, -58};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {57, 58, -15, -14};
//+
Plane Surface(20) = {20};
//+
Physical Surface(2) = {20, 18, 16, 14, 12, 10, 7, 8, 4, 2, 1, 3, 5, 6, 9, 11, 13, 15, 17, 19};
//+
Physical Curve(0) = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
//+
Physical Curve(11) = {15, 16};
//+
Physical Curve(13) = {1, 2};

// Meshing

Mesh.Algorithm = 3;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 1.;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
