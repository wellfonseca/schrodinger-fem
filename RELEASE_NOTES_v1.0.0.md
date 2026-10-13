# schrodinger-fem v1.0.0

First public release. A finite element solver for the stationary Schrödinger
eigenvalue problem in one, two and three spatial dimensions, as a module of the
[Elmer](https://www.csc.fi/web/elmer) package, with an independent Python
implementation for code-to-code verification.

## What is included

* **`src/`** — the Fortran solver module `SchrodingerSolver.F90` and
  `elmer_schrodinger_files.def`, an Apptainer recipe that builds Elmer with the
  module already in place for machines with no Elmer installation.
* **`cases/`** — seven self-contained Elmer benchmark cases, each with its Gmsh
  mesh source, mesh, `.sif` input file, the exact commands and the output of the
  run performed here.
* **`python/`** — the independent implementation, with a dense LAPACK route and
  a sparse shift-and-invert ARPACK route, plus matrix export for PETSc/SLEPc.
* **`docs/`** — model and solver manuals in English and Portuguese.
* **`verification/`** — the scripts that reproduce every numerical table of the
  accompanying manuscript.

## The seven benchmark cases

| case | dimension | domain | potential | elements | $n_f$ |
|---|---|---|---|---|---|
| `1D_poco_infinito` | 1D | $[-20,20]$ | $V=0$ | 250 | 249 |
| `1D_harmonico` | 1D | $[-10,10]$ | $V=x^2$ | 400 | 399 |
| `1D_morse` | 1D | $[0,10]$ | $V=10(e^{-2a(x-x_0)}-2e^{-a(x-x_0)})$ | 400 | 399 |
| `1D_poco_duplo` | 1D | $[-10,10]$ | $V=0.01(x^2-4)^2$ | 400 | 399 |
| `1D_senoidal` | 1D | $[-10,10]$ | $V=0.5\cos(2\pi x/10)$ | 400 | 399 |
| `2D_harmonico` | 2D | $[-10,10]^2$ | $V=0.25(x^2+y^2)$ | $40\times40$ | 1521 |
| `3D_harmonico` | 3D | $[-6,6]^3$ | $V=0.25(x^2+y^2+z^2)$ | $14^3$ | 2197 |

## Verification

Four independent checks, all reproducible from this repository:

1. **Analytic benchmarks** — the seven cases above against closed-form spectra.
2. **Convergence order** — $\mathcal{O}(h^2)$ for eigenvalues, in 1D, 2D and 3D.
3. **Code-to-code** — the Elmer module and the Python code agree to about
   $4\times10^{-13}$ in relative terms on identical meshes.
4. **Structural invariants** — degeneracies, mass-orthonormality, and the virial
   relation.

## Documented failure modes

Two problems that are easy to encounter and hard to diagnose are documented in
full, with numerical evidence:

* **The diagonal-shift treatment of Dirichlet data is invalid for an
  eigenproblem.** It introduces spurious eigenvalues, at exactly $E=1$ if both
  diagonals are set to one, or at a distorted value if only the stiffness
  diagonal is modified. The module eliminates the constrained degrees of
  freedom instead, which is the correct procedure. An earlier version of the
  accompanying analysis argued for this via a block-matrix argument that happens
  to hold for exponentially decaying eigenfunctions and fails for the infinite
  well, where the residual of the discarded rows reaches 8 % of the amplitude.
* **The module assembles dense local matrices**, so `ElmerSolver_mpi` returns
  wrong results with more than one rank. Use the serial `ElmerSolver`; see
  `docs/PARALELISMO.md`.

## Scope and limitations

* Arbitrary scalar potentials in one dimension; **separable** potentials in two
  and three dimensions.
* Homogeneous Dirichlet boundary conditions only.
* No distributed timings are claimed: the exported matrices can be solved with
  PETSc/SLEPc, but a parallel sparse direct factorization was not available in
  the installation used.

## License

GPL-3.0-or-later, matching Elmer. See `LICENSE` and `NOTICE`.
