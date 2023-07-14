// Constants //
// Mesh.SubdivisionAlgorithm = 1;
length = 12*0.0254;
height = 0.0262; // 0.0265
lc1 = 0.3e-3; // 2e-3
lc2 = 2.0e-5; // Normally lc1/75
critical1 = 0.005*length;
critical2 = 0.005*length; // 0.03*length
nProx = 3 * critical1;
// Right Points //
Point(7) = {length/2, 0, 0, lc1};
Point(8) = {nProx, 0, 0, lc1};
Point(9) = {critical1, 0, 0, lc2};
Point(10) = {critical1, height, 0, lc2};
Point(11) = {nProx, height, 0, lc1};
Point(12) = {length/2, height, 0, lc1};
// Center Points //
Point(13) = {0, 0, 0, lc2};
Point(14) = {0, height, 0, lc2};
Point(15) = {critical1, height/2, 0, lc1};
Point(16) = {0, height/2, 0, lc1};
// Right Lines //
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 7};
// Center Lines //
Line(13) = {13, 16};
Line(14) = {16, 14};
Line(15) = {14, 10};

Line(17) = {9, 15};
Line(18) = {15, 10};
Line(19) = {9, 13};



// Line Loop //
Line Loop(2) = {7, 8, 17, 18, 10, 11, 12}; // Right Plane
Line Loop(3) = {13, 14, 15, -18, -17, 19}; // Center Plane
// Surface //
Plane Surface(2) = {2};
Plane Surface(3) = {3};
