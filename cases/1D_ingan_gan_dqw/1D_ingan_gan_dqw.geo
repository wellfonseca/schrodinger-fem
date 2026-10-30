// ---------------------------------------------------------------------------
//  1D_ingan_gan_dqw  --  symmetric InGaN/GaN double quantum well
//
//  Case reproduced from
//    En-nadir, E. et al., J. Comput. Electron. (2022),
//    doi:10.1007/s10825-021-01830-4
//
//  Structure: well(3 nm) | barrier(3 nm) | well(3 nm), total 9 nm,
//  embedded in a 16 nm domain so that the Dirichlet walls do not touch it.
//
//     gmsh 1D_ingan_gan_dqw.geo -2 -o 1D_ingan_gan_dqw.msh
//     ElmerGrid 14 2 1D_ingan_gan_dqw.msh -out mesh
// ---------------------------------------------------------------------------
L  = 16.0;
n  = 640;
x0 = 0.0;
x1 = x0 + L;

Point(1) = {x0, 0, 0, 1.0};
Point(2) = {x1, 0, 0, 1.0};
Line(1)  = {1, 2};
Transfinite Line {1} = n + 1;
Physical Point(1) = {1, 2};
Physical Line(1)  = {1};
