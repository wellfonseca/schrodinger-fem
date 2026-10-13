// ---------------------------------------------------------------------------
//  1D_infinite_well  --  structured 1D mesh (GMSH)
//  L = 40.0  n = 250  origin x0 = -20.0
//
//  Generate the .msh and convert it to the Elmer format:
//     gmsh 1D_infinite_well.geo -2 -o 1D_infinite_well.msh
//     ElmerGrid 14 2 1D_infinite_well.msh -out mesh
//  ElmerGrid reads the Gmsh .msh with code 14 and writes the ElmerSolver
//  format (code 2) into the mesh/ directory.
// ---------------------------------------------------------------------------
L  = 40.0;
n  = 250;
x0 = -20.0;
x1 = x0 + L;

Point(1) = {x0, 0, 0, 1.0};
Point(2) = {x1, 0, 0, 1.0};
Line(1)  = {1, 2};
Transfinite Line {1} = n + 1;
Physical Point(1) = {1, 2};
Physical Line(1)  = {1};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm = 1;
