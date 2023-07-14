// Constants //
// Mesh.SubdivisionAlgorithm = 1;
Mesh.Algorithm = 5;
length = 12*0.0254;
height = 0.0262; // 0.0265
lc1 = 300e-6; // 2e-3
lc2 = 16.67e-6; // Normally lc1/75
critical1 = 0.005*length;
critical2 = 0.005*length; // 0.03*length
nProx = 3 * critical1;
// Left Points //
Point(1) = {-length/2, 0, 0, lc1};
Point(2) = {-nProx, 0, 0, lc1};
Point(3) = {-critical1, 0, 0, lc2};
Point(4) = {-critical1, height, 0, lc2};
Point(5) = {-nProx, height, 0, lc1};
Point(6) = {-length/2, height, 0, lc1};
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
Point(16) = {-critical1, height/2, 0, lc1};

// Left Lines //
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
// Right Lines //
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 7};
// Center Lines //
Line(13) = {3, 13};
Line(14) = {13, 9};
Line(15) = {10, 14};
Line(16) = {14, 4};

Line(17) = {9, 15};
Line(18) = {15, 10};
Line(19) = {4, 16};
Line(20) = {16, 3};


// Line Loop //
Line Loop(1) = {1, 2, -20, -19, 4, 5, 6}; // Left Plane
Line Loop(2) = {7, 8, 17, 18, 10, 11, 12}; // Right Plane
Line Loop(3) = {13, 14, 17, 18, 15, 16, 19, 20}; // Center Plane
// Surface //
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
