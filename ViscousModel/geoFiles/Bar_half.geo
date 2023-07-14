// Constants //
// Mesh.SubdivisionAlgorithm = 1;
length = 12*0.0254;
height = 0.0265;
lc1 = 0.002;
lc2 = lc1/75; // Normally lc1/75
critical1 = 0.05*length;
critical2 = 0.05*length; // 0.03*length
// General Points //
Point(4) = {0, 0, 0, lc2};
Point(5) = {critical1, 0, 0, lc2};
Point(6) = {length/5, 0, 0, lc1};
Point(7) = {length/2, 0, 0, lc1};
Point(8) = {length/2, height, 0, lc1};
Point(9) = {length/5, height, 0, lc1};
Point(10) = {critical2, height, 0, lc2};
Point(11) = {0, height, 0, lc2};
// Lines //
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 4};
// Line Loop //
Line Loop(1) = {4, 5, 6, 7, 8, 9, 10, 11};
// Surface //
Plane Surface(1) = {1};
