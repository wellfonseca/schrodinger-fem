#!/usr/bin/env python3
"""
Generate the Gmsh .geo sources for the benchmark cases.

Structured meshes (transfinite + recombine) matching what ElmerSolver expects
in 1D, 2D and 3D. The .msh is written in version 2.2, which ElmerGrid reads
with the format code 14.
"""
import os

# case: (dim, L, n, x0)
CASES = {
    '1D_infinite_well': (1, 40.0, 250, -20.0),
    '1D_harmonic':      (1, 20.0, 400, -10.0),
    '1D_morse':         (1, 10.0, 400,   0.0),
    '1D_double_well':   (1, 20.0, 400, -10.0),
    '1D_sinusoidal':    (1, 20.0, 400, -10.0),
    '2D_harmonic':      (2, 20.0,  40, -10.0),
    '3D_harmonic':      (3, 12.0,  14,  -6.0),
}

GEO_1D = """// ---------------------------------------------------------------------------
//  {case}  --  structured 1D mesh (GMSH)
//  L = {L}  n = {n}  origin x0 = {x0}
//
//  Generate the .msh and convert it to the Elmer format:
//     gmsh {case}.geo -2 -o {case}.msh
//     ElmerGrid 14 2 {case}.msh -out mesh
//  ElmerGrid reads the Gmsh .msh with code 14 and writes the ElmerSolver
//  format (code 2) into the mesh/ directory.
// ---------------------------------------------------------------------------
L  = {L};
n  = {n};
x0 = {x0};
x1 = x0 + L;

Point(1) = {{x0, 0, 0, 1.0}};
Point(2) = {{x1, 0, 0, 1.0}};
Line(1)  = {{1, 2}};
Transfinite Line {{1}} = n + 1;
Physical Point(1) = {{1, 2}};
Physical Line(1)  = {{1}};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm = 1;
"""

GEO_2D = """// ---------------------------------------------------------------------------
//  {case}  --  structured 2D quadrilateral mesh (GMSH)
//  L = {L}  n = {n} per direction  origin ({x0},{x0})
//
//  Generate the .msh and convert it:
//     gmsh {case}.geo -2 -o {case}.msh
//     ElmerGrid 14 2 {case}.msh -out mesh
// ---------------------------------------------------------------------------
L  = {L};
n  = {n};
x0 = {x0};
x1 = x0 + L;

Point(1) = {{x0, x0, 0, 1.0}};
Point(2) = {{x1, x0, 0, 1.0}};
Point(3) = {{x1, x1, 0, 1.0}};
Point(4) = {{x0, x1, 0, 1.0}};

Line(1) = {{1, 2}};
Line(2) = {{2, 3}};
Line(3) = {{3, 4}};
Line(4) = {{4, 1}};

Line Loop(1) = {{1, 2, 3, 4}};
Plane Surface(1) = {{1}};

Transfinite Line {{1, 2, 3, 4}} = n + 1;
Transfinite Surface {{1}} = {{1, 2, 3, 4}};
Recombine Surface {{1}};

Physical Surface(1) = {{1}};
Physical Line(1) = {{1, 2, 3, 4}};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm = 1;
Mesh.RecombineAll = 1;
"""

GEO_3D = """// ---------------------------------------------------------------------------
//  {case}  --  structured 3D hexahedral mesh (GMSH)
//  L = {L}  n = {n} per direction  origin ({x0},{x0},{x0})
//
//  Generate the .msh and convert it:
//     gmsh {case}.geo -2 -o {case}.msh
//     ElmerGrid 14 2 {case}.msh -out mesh
// ---------------------------------------------------------------------------
L  = {L};
n  = {n};
h  = L/2;
x0 = {x0};
x1 = x0 + L;

Point(1) = {{x0, x0, x0, 1.0}};
Point(2) = {{x1, x0, x0, 1.0}};
Point(3) = {{x1, x1, x0, 1.0}};
Point(4) = {{x0, x1, x0, 1.0}};
Point(5) = {{x0, x0, x1, 1.0}};
Point(6) = {{x1, x0, x1, 1.0}};
Point(7) = {{x1, x1, x1, 1.0}};
Point(8) = {{x0, x1, x1, 1.0}};

Line(1) = {{1, 2}};  Line(2) = {{2, 3}};  Line(3) = {{3, 4}};  Line(4) = {{4, 1}};
Line(5) = {{5, 6}};  Line(6) = {{6, 7}};  Line(7) = {{7, 8}};  Line(8) = {{8, 5}};
Line(9) = {{1, 5}};  Line(10) = {{2, 6}}; Line(11) = {{3, 7}}; Line(12) = {{4, 8}};

Line Loop(1) = {{1, 2, 3, 4}};   Plane Surface(1) = {{1}};
Line Loop(2) = {{5, 6, 7, 8}};   Plane Surface(2) = {{2}};
Line Loop(3) = {{1, 10, -5, -9}};  Plane Surface(3) = {{3}};
Line Loop(4) = {{2, 11, -6, -10}}; Plane Surface(4) = {{4}};
Line Loop(5) = {{3, 12, -7, -11}}; Plane Surface(5) = {{5}};
Line Loop(6) = {{4, 9, -8, -12}};  Plane Surface(6) = {{6}};

Surface Loop(1) = {{1, 2, 3, 4, 5, 6}};
Volume(1) = {{1}};

Transfinite Line {{1,2,3,4,5,6,7,8,9,10,11,12}} = n + 1;
Transfinite Surface {{1}} = {{1, 2, 3, 4}};
Transfinite Surface {{2}} = {{5, 6, 7, 8}};
Transfinite Surface {{3}} = {{1, 2, 6, 5}};
Transfinite Surface {{4}} = {{2, 3, 7, 6}};
Transfinite Surface {{5}} = {{3, 4, 8, 7}};
Transfinite Surface {{6}} = {{4, 1, 5, 8}};
Recombine Surface {{1, 2, 3, 4, 5, 6}};
Transfinite Volume {{1}} = {{1, 2, 3, 4, 5, 6, 7, 8}};
Recombine Volume {{1}};

Physical Volume(1)  = {{1}};
Physical Surface(1) = {{1, 2, 3, 4, 5, 6}};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm = 1;
Mesh.RecombineAll = 1;
"""

TEMPLATE = {1: GEO_1D, 2: GEO_2D, 3: GEO_3D}

if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    for case, (dim, L, n, x0) in CASES.items():
        txt = TEMPLATE[dim].format(case=case, L=L, n=n, x0=x0)
        out = os.path.join(here, case, '%s.geo' % case)
        with open(out, 'w', encoding='utf-8') as f:
            f.write(txt)
        print("  %-20s %s.geo  (%dD, L=%g, n=%d)" % (case, case, dim, L, n))
