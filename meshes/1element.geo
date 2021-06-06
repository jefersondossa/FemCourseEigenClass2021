// Gmsh project created on Mon May 24 18:40:06 2021
//+
Point(1) = {0, 0, 0, 0.1};
Point(2) = {1, 0, 0, 0.1};

Line(1) = {1, 2};

//Transfinite Curve {1} = 11 Using Progression 1;

Physical Curve("1") = {1};//Bottom Line

