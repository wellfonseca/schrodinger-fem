// ---------------------------------------------------------------------------
//  2D_harmonico  --  malha 2D estruturada de quadrilateros (GMSH)
//  L = 20.0  n = 40 por direcao   origem (-10.0,-10.0)
//
//  Gerar o .msh e converter:
//     gmsh 2D_harmonico.geo -2 -o 2D_harmonico.msh
//     ElmerGrid 14 2 2D_harmonico.msh -out mesh
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
