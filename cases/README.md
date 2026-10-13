# Benchmark cases

Seven self-contained Elmer cases. Each folder holds the mesh (`.geo` Gmsh
source and `.msh` mesh), the Elmer input file `caso.sif`, the exact commands in
`command.txt`, and the output of the run performed here in `run.log`.

## Running a case

```bash
cd 1D_harmonico
ElmerGrid 14 2 1D_harmonico.msh -out mesh    # 14 = Gmsh input, 2 = Elmer output
ElmerSolver caso.sif
```

The eigenvalues are reported in `run.log` and a `.vtu` file is written to
`mesh/` for ParaView. With the container built from `src/`:

```bash
apptainer exec ../elmer_schrodinger.sif ElmerSolver caso.sif
```

## The cases

| case | dimension | domain | potential | elements | $n_f$ |
|---|---|---|---|---|---|
| `1D_poco_infinito` | 1D | $[-20,20]$ | $V=0$ | 250 | 249 |
| `1D_harmonico` | 1D | $[-10,10]$ | $V=x^2$ | 400 | 399 |
| `1D_morse` | 1D | $[0,10]$ | $V=10(e^{-2a(x-x_0)}-2e^{-a(x-x_0)})$ | 400 | 399 |
| `1D_poco_duplo` | 1D | $[-10,10]$ | $V=0.01(x^2-4)^2$ | 400 | 399 |
| `1D_senoidal` | 1D | $[-10,10]$ | $V=0.5\cos(2\pi x/10)$ | 400 | 399 |
| `2D_harmonico` | 2D | $[-10,10]^2$ | $V=0.25(x^2+y^2)$ | $40\times40$ | 1521 |
| `3D_harmonico` | 3D | $[-6,6]^3$ | $V=0.25(x^2+y^2+z^2)$ | $14^3$ | 2197 |

$n_f$ is the number of free degrees of freedom after the Dirichlet elimination.

## Regenerating the meshes

`gerar_malha.py` in each folder regenerates the `.geo` and `.msh` with Gmsh;
`../../cases/gera_todas_malhas.sh` does all seven at once.

## Note on distributed runs

Use the serial `ElmerSolver`. `ElmerSolver_mpi` with more than one rank returns
wrong results, because the module assembles dense local matrices. See
`../docs/PARALELISMO.md`.
