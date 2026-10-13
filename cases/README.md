# Benchmark cases

Seven self-contained Elmer cases. Each folder holds only the case itself:

| file | purpose |
|---|---|
| `<case>.geo` | Gmsh geometry source, structured transfinite mesh |
| `<case>.msh` | the mesh in Gmsh format, as committed |
| `case.sif` | the Elmer input file |
| `command.txt` | the exact commands to run the case |

The mesh is generated with Gmsh from the `.geo` source; both the source and the
resulting `.msh` are committed, so nothing has to be regenerated to run a case.

## Running a case

```bash
cd 1D_harmonic
ElmerGrid 14 2 1D_harmonic.msh -out mesh    # 14 = Gmsh input, 2 = Elmer output
ElmerSolver case.sif
```

The eigenvalues are printed to standard output and a `.vtu` file is written to
`mesh/` for ParaView. With the container built from `../src/`:

```bash
apptainer exec elmer_schrodinger.sif ElmerSolver case.sif
```

## Regenerating the meshes

The meshes come from Gmsh; the `.geo` sources are the input and the `.msh`
files are the result:

```bash
cd 1D_harmonic
gmsh 1D_harmonic.geo -2 -o 1D_harmonic.msh
ElmerGrid 14 2 1D_harmonic.msh -out mesh
```

The `.geo` files use transfinite lines and recombined surfaces so that the
structured mesh matches what ElmerSolver expects, and `Mesh.MshFileVersion = 2.2`
so that ElmerGrid reads the result with format code 14.

## The cases

| case | dimension | domain | potential | elements | $n_f$ |
|---|---|---|---|---|---|
| `1D_infinite_well` | 1D | $[-20,20]$ | $V=0$ | 250 | 249 |
| `1D_harmonic` | 1D | $[-10,10]$ | $V=x^2$ | 400 | 399 |
| `1D_morse` | 1D | $[0,10]$ | $V=10(e^{-2a(x-x_0)}-2e^{-a(x-x_0)})$, $a=0.5$, $x_0=3$ | 400 | 399 |
| `1D_double_well` | 1D | $[-10,10]$ | $V=0.01(x^2-4)^2$ | 400 | 399 |
| `1D_sinusoidal` | 1D | $[-10,10]$ | $V=0.5\cos(2\pi x/10)$ | 400 | 399 |
| `2D_harmonic` | 2D | $[-10,10]^2$ | $V=0.25(x^2+y^2)$ | $40\times40$ | 1521 |
| `3D_harmonic` | 3D | $[-6,6]^3$ | $V=0.25(x^2+y^2+z^2)$ | $14^3$ | 2197 |

$n_f$ is the number of free degrees of freedom after the Dirichlet elimination.

## Note on distributed runs

Use the serial `ElmerSolver`. `ElmerSolver_mpi` with more than one rank returns
wrong results, because the module assembles dense local matrices. See
`../docs/PARALLELISM.md`.
