// Constants //
length = 12*0.0254;
height = 0.0265;
lc1 = 0.002;
lc2 = lc1/75; // Normally lc1/75
critical1 = 0.03*length;
critical2 = 0.03*length; // 0.03*length
// General Points //
Point(1) = {-length/2, 0, 0, lc1};
Point(2) = {-length/5, 0, 0, lc1};
Point(3) = {-critical1, 0, 0, lc2};
Point(4) = {0, 0, 0, lc2};
Point(5) = {critical1, 0, 0, lc2};
Point(6) = {length/5, 0, 0, lc1};
Point(7) = {length/2, 0, 0, lc1};
Point(8) = {length/2, height, 0, lc1};
Point(9) = {length/5, height, 0, lc1};
Point(10) = {critical2, height, 0, lc2};
Point(11) = {0, height, 0, lc2};
Point(12) = {-critical2, height, 0, lc2};
Point(13) = {-length/5, height, 0, lc1};
Point(14) = {-length/2, height, 0, lc1};
// Lines //
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 1};
// Line Loop //
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
// Surface //
Plane Surface(1) = {1};
