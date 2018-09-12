lc = 1;

L = 2;
H = 1;
W = 0.01;

Nx = 100;
Ny = 50;
Nz = 1;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, H, 0, lc};
Point(4) = {0, H, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line {1, 3} = Nx Using Progression 1;
Transfinite Line {2, 4} = Ny Using Progression 1;

Line Loop(11) = {1, 2, 3, 4};
Plane Surface(21) = {11};
Transfinite Surface {21};
Recombine Surface {21};

Extrude {0, 0, W} {
    Surface{21}; Layers{Nz}; Recombine;
}

Physical Surface("inlet") = {42};
Physical Surface("outlet") = {34};
Physical Surface("topAndBottom") = {38, 30};
Physical Surface("frontAndBack") = {43, 21};
Physical Volume("cavity") = {1};

