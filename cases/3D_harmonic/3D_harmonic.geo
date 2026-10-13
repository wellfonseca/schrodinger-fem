// ---------------------------------------------------------------------------
//  3D_harmonic  --  structured 3D hexahedral mesh (GMSH)
//  L = 12.0  n = 14 per direction  origin (-6.0,-6.0,-6.0)
//
//  Generate the .msh and convert it:
//     gmsh 3D_harmonic.geo -2 -o 3D_harmonic.msh
//     ElmerGrid 14 2 3D_harmonic.msh -out mesh
// ---------------------------------------------------------------------------
L  = 12.0;
n  = 14;
h  = L/2;
x0 = -6.0;
x1 = x0 + L;

Point(1) = {x0, x0, x0, 1.0};
Point(2) = {x1, x0, x0, 1.0};
Point(3) = {x1, x1, x0, 1.0};
Point(4) = {x0, x1, x0, 1.0};
Point(5) = {x0, x0, x1, 1.0};
Point(6) = {x1, x0, x1, 1.0};
Point(7) = {x1, x1, x1, 1.0};
Point(8) = {x0, x1, x1, 1.0};

Line(1) = {1, 2};  Line(2) = {2, 3};  Line(3) = {3, 4};  Line(4) = {4, 1};
Line(5) = {5, 6};  Line(6) = {6, 7};  Line(7) = {7, 8};  Line(8) = {8, 5};
Line(9) = {1, 5};  Line(10) = {2, 6}; Line(11) = {3, 7}; Line(12) = {4, 8};

Line Loop(1) = {1, 2, 3, 4};   Plane Surface(1) = {1};
Line Loop(2) = {5, 6, 7, 8};   Plane Surface(2) = {2};
Line Loop(3) = {1, 10, -5, -9};  Plane Surface(3) = {3};
Line Loop(4) = {2, 11, -6, -10}; Plane Surface(4) = {4};
Line Loop(5) = {3, 12, -7, -11}; Plane Surface(5) = {5};
Line Loop(6) = {4, 9, -8, -12};  Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12} = n + 1;
Transfinite Surface {1} = {1, 2, 3, 4};
Transfinite Surface {2} = {5, 6, 7, 8};
Transfinite Surface {3} = {1, 2, 6, 5};
Transfinite Surface {4} = {2, 3, 7, 6};
Transfinite Surface {5} = {3, 4, 8, 7};
Transfinite Surface {6} = {4, 1, 5, 8};
Recombine Surface {1, 2, 3, 4, 5, 6};
Transfinite Volume {1} = {1, 2, 3, 4, 5, 6, 7, 8};
Recombine Volume {1};

Physical Volume(1)  = {1};
Physical Surface(1) = {1, 2, 3, 4, 5, 6};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm = 1;
Mesh.RecombineAll = 1;
