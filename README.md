# schrodinger-fem

A finite element solver for the **stationary Schrödinger eigenvalue problem** in
one, two and three spatial dimensions, implemented as a module of the
[Elmer](https://www.csc.fi/web/elmer) finite element package, together with an
independent Python implementation used for code-to-code verification.

The problem solved is

$$-\nabla^2\psi + V(\mathbf{r})\,\psi = E\,\psi \quad \text{in } \Omega, \qquad \psi = 0 \quad \text{on } \partial\Omega,$$

in the unit convention $\hbar = m = 1$. The weak form leads to the generalized
symmetric eigenproblem

$$A\,\boldsymbol{\psi} = E\,M\,\boldsymbol{\psi}, \qquad A = K + M_V,$$

with $K_{ij} = \int \nabla N_i \cdot \nabla N_j$, $M_{ij} = \int N_i N_j$ and
$(M_V)_{ij} = \int N_i V N_j$, assembled with multilinear Lagrange elements.

## What this package contains

| directory | contents |
|---|---|
| `src/` | the Elmer solver module `SchrodingerSolver.F90` and the build recipe |
| `cases/` | seven Elmer benchmark cases, each self-contained and runnable |
| `python/` | the independent Python implementation and its cases |
| `article/` | the LaTeX source and the compiled manuscript |
| `docs/` | model manual and solver manual, in English and Portuguese |
| `verification/` | scripts that reproduce the tables and figures of the article |

## Quick start — the Elmer module

The solver is a standard Elmer solver module: it is compiled against an existing
Elmer installation and dropped next to the other solvers. Full instructions are
in [`docs/Manual_Solver_Schrodinger_EN.pdf`](docs/Manual_Solver_Schrodinger_EN.pdf);
in outline, from a built Elmer tree:

```bash
cd $ELMER_HOME/elmersolver
cp /path/to/src/SchrodingerSolver.F90 .
# add SchrodingerSolver to the module list of the elmersolver Makefile
make
```

The recipe in `src/elmer_schrodinger_files.def` builds an Apptainer image with
Elmer and the module already in place, for a machine with no Elmer installation.

## Quick start — a benchmark case

Each case under `cases/` is self-contained. For example:

```bash
cd cases/1D_harmonico
ElmerGrid 14 2 1D_harmonico.msh -out mesh   # Gmsh -> Elmer
ElmerSolver caso.sif
```

The computed eigenvalues are written to `run.log` and a `.vtu` file is left in
`mesh/` for ParaView. `command.txt` in each case folder records the exact
commands, and `run.log` records the result obtained here.

The seven cases are:

| case | domain | potential | elements |
|---|---|---|---|
| `1D_poco_infinito` | $[-20,20]$ | $V=0$ | 250 |
| `1D_harmonico` | $[-10,10]$ | $V=x^2$ | 400 |
| `1D_morse` | $[0,10]$ | $V=10(e^{-2a(x-x_0)}-2e^{-a(x-x_0)})$, $a=0.5$, $x_0=3$ | 400 |
| `1D_poco_duplo` | $[-10,10]$ | $V=0.01(x^2-4)^2$ | 400 |
| `1D_senoidal` | $[-10,10]$ | $V=0.5\cos(2\pi x/10)$ | 400 |
| `2D_harmonico` | $[-10,10]^2$ | $V=0.25(x^2+y^2)$ | $40\times40$ |
| `3D_harmonico` | $[-6,6]^3$ | $V=0.25(x^2+y^2+z^2)$ | $14^3$ |

## Quick start — the Python implementation

```bash
cd python
python3 fem_schrodinger.py --case casos/1D_harmonico/run.py
```

`fem_schrodinger.py` assembles the matrices and solves the dense generalized
eigenproblem with LAPACK; `eigsolver_esparso.py` solves the same problem by
sparse shift-and-invert Arnoldi iteration (ARPACK), which extends the accessible
problem size by three orders of magnitude in the number of degrees of freedom.

Requires Python 3 with `numpy` and `scipy`.

## Units and the effective-mass interpretation

The solver uses a unit coefficient on the Laplacian. For a semiconductor
heterostructure with effective mass $m^*$ and characteristic length $L$, the
physical energy unit is

$$\mathcal{E} = \frac{\hbar^2}{2m^*L^2}.$$

For GaAs ($m^* \approx 0.067\,m_e$) and $L = 10$ nm this gives
$\mathcal{E} \approx 5.7$ meV. Multiplying a computed eigenvalue by
$\mathcal{E}$ gives the physical energy; see
[`docs/Manual_Model_Schrodinger_EN.pdf`](docs/Manual_Model_Schrodinger_EN.pdf).

## Scope and limitations

* **Potentials.** Arbitrary scalar potentials $V(x)$ in one dimension.
  In two and three dimensions the implementation supports **separable**
  potentials, $V = \sum_k V_k(x_k)$. A general non-separable $V(\mathbf{r})$ is
  not supported by the current assembly.
* **Boundary conditions.** Homogeneous Dirichlet only. The diagonal-shift
  treatment familiar from source problems is **invalid here**: it introduces
  spurious eigenvalues. The code eliminates the constrained degrees of freedom,
  which is the correct procedure. See Section 3.3 of the manuscript.
* **Parallel execution.** The assembly is written for a partitioned mesh, but
  the module currently assembles dense local matrices, so `ElmerSolver_mpi`
  returns wrong results with more than one rank. This is documented in
  `docs/PARALELISMO.md`; use the serial `ElmerSolver` for correct results.
* **The distributed route.** The exported matrices can be solved with
  PETSc/SLEPc, but a parallel sparse direct factorization was not available in
  the installation used, so no distributed timings are claimed.

## Verification

The implementation is verified in four independent ways, all reproducible from
this repository:

1. **Analytic benchmarks** — the seven cases above against closed-form spectra
   (`verification/confere.py`).
2. **Convergence order** — $\mathcal{O}(h^2)$ for eigenvalues, verified in 1D, 2D
   and 3D (`verification/conv1d.py`, `conv2d_int.py`, `conv3d_int.py`).
3. **Code-to-code** — the Elmer module and the Python code agree to
   $\sim4\times10^{-13}$ relative on identical meshes.
4. **Structural invariants** — degeneracies, mass-orthonormality, and the virial
   relation (`verification/orto.py`, `verification/virial4.py`).

## Citing

If you use this software, please cite the manuscript and the archived release.
See [`CITATION.cff`](CITATION.cff).

## License

GPL-3.0-or-later, matching Elmer. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
