// Parameters
L = 0.5e-3; // size of the square (1 mm)
length_scale = 4.e-6;
factor = 2;
//refinement_size = slit_size; // Local refinement size (0.001 mm)
//refinement_length = 0.5; // Length of the refinement region across the middle

//+
//SetFactory("OpenCASCADE");
//+
//Rectangle(1) = {-L, -L, 0, 2*L, 2*L, 0};
//+
//Rectangle(2) = {-L-.1*L, 0, 0, L+.1*L, length_scale/4., 0};
//+
//BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }

Point(1) = {-L, -L, 0., L/4.};
Point(2) = {L, -L, 0., L/4.};
Point(3) = {L, 0., 0., length_scale/factor};
Point(4) = {0., 0., 0., length_scale/factor};
Point(5) = {-L, 0., 0., L/4.};

Point(6) = {-L, length_scale/4., 0., L/4.};
Point(7) = {0., length_scale/4., 0., length_scale/factor};
Point(8) = {L, length_scale/4., 0., length_scale/factor};
Point(9) = {L, L, 0., L/4.};
Point(10) = {-L, L, 0., L/4.};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 1};

//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 6};

//+
Line(11) = {3, 8};
//+
Line(12) = {7, 4};

//+
Curve Loop(1) = {1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};

//+
Curve Loop(2) = {6, 7, 8, 9, 10};
//+
Plane Surface(2) = {2};

//+
Curve Loop(3) = {3, -12, 7, -11};
//+
Plane Surface(3) = {3};

Physical Surface(2) = {1, 3, 2};

Physical Curve(13) = {1};

Physical Curve(11) = {9};

Mesh.Algorithm = 3;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 1.;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";

