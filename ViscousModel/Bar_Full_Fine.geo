// Constants //
// Mesh.SubdivisionAlgorithm = 1;
Mesh.Algorithm = 5; // Was 6
length = 12*0.0254; // 12*0.0254
height = 0.0210; // 0.0265
lc1 = 300e-6; // 300e-6
lc2 = 16.67e-6; // Normally 16.67e-6
halfModifier = 1; // Simply mesh if = 2, normal if = 1
centerLower = (1/4) / halfModifier; // Normally 1/4, must be lower than centerHigher
centerHigher = (1/3) / halfModifier; // Normally 1/3, must be higher than centerLower 
critical1 = 0.005*length; // 0.005*length
critical2 = 0.005*length; // 0.005*length
nProx = 3 * critical1; // 3 * critical1

// Left Points //
Point(1) = {-length/2, 0, 0, lc1};
Point(2) = {-nProx, 0, 0, lc1};
Point(3) = {-critical1, 0, 0, lc2};
Point(4) = {-critical1, height * centerLower, 0, lc2};
Point(5) = {-critical1, height * centerHigher, 0, lc1};
Point(6) = {-critical1, height, 0, lc2};
Point(7) = {-nProx, height, 0, lc1};
Point(8) = {-length/2, height, 0, lc1};

// Right Points //
Point(9) = {length/2, 0, 0, lc1};
Point(10) = {nProx, 0, 0, lc1};
Point(11) = {critical1, 0, 0, lc2};
Point(12) = {critical1, height * centerLower, 0, lc2};
Point(13)  ={critical1, height * centerHigher, 0, lc1};
Point(14) = {critical1, height, 0, lc2};
Point(15) = {nProx, height, 0, lc1};
Point(16) = {length/2, height, 0, lc1};

// Center Points //
Point(17) = {0, 0, 0, lc2};
Point(18) = {0, height, 0, lc2};

// Left Lines //
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

// Right Lines //
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 9};

// Center Lines //
Line(17) = {3, 17};
Line(18) = {17, 11};
Line(19) = {14, 18};
Line(20) = {18, 6};

// Line Loop //
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8}; // Left Plane
Line Loop(2) = {9, 10, 11, 12, 13, 14, 15, 16}; // Right Plane
Line Loop(3) = {17, 18, 11, 12, 13, 19, 20, -5, -4, -3}; // Center Plane

// Surface //
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
