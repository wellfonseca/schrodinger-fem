// ---------------------------------------------------------------------------
//  1D_poco_infinito  --  malha 1D estruturada (GMSH)
//  L = 40.0  n = 250  origem x0 = -20.0
//
//  Gerar o .msh e converter para o formato Elmer:
//     gmsh 1D_poco_infinito.geo -2 -o 1D_poco_infinito.msh
//     ElmerGrid 14 2 1D_poco_infinito.msh -out mesh
//  O ElmerGrid le o .msh (Gmsh) com o codigo 14 e escreve o formato
//  ElmerSolver (codigo 2) no diretorio 1D_poco_infinito/.
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
