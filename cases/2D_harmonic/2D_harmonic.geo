// ---------------------------------------------------------------------------
//  2D_harmonic  --  structured 2D quadrilateral mesh (GMSH)
//  L = 20.0  n = 40 per direction  origin (-10.0,-10.0)
//
//  Generate the .msh and convert it:
//     gmsh 2D_harmonic.geo -2 -o 2D_harmonic.msh
//     ElmerGrid 14 2 2D_harmonic.msh -out mesh
// ---------------------------------------------------------------------------
L  = 20.0;
n  = 40;
x0 = -10.0;
x1 = x0 + L;

Point(1) = {x0, x0, 0, 1.0};
Point(2) = {x1, x0, 0, 1.0};
Point(3) = {x1, x1, 0, 1.0};
Point(4) = {x0, x1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 2, 3, 4} = n + 1;
Transfinite Surface {1} = {1, 2, 3, 4};
Recombine Surface {1};

Physical Surface(1) = {1};
Physical Line(1) = {1, 2, 3, 4};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm = 1;
Mesh.RecombineAll = 1;
