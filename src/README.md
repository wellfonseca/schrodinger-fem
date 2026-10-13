# src — the Elmer solver module

`SchrodingerSolver.F90` is a standard Elmer solver module. It reads a potential
from the `.sif` file, assembles the stiffness, mass and potential-weighted mass
matrices, eliminates the constrained degrees of freedom and solves the resulting
generalized eigenproblem.

## Adding it to an Elmer installation

The module must be compiled together with the other Elmer solvers, because it
uses their module interfaces:

```bash
cd /path/to/elmerfem/elmersolver
cp /path/to/src/SchrodingerSolver.F90 .
```

Then add `SchrodingerSolver` to the list of solver objects in the
`elmersolver/Makefile` and rebuild:

```bash
make
make install
```

The solver then appears among the native solvers and can be selected in a
`.sif` file with `Solver = "SchrodingerSolver"`.

## Building without an Elmer installation

`elmer_schrodinger_files.def` is an Apptainer definition that builds an image
containing Elmer and the module already in place, for a machine with no Elmer
installation and no root access:

```bash
apptainer build elmer_schrodinger.sif elmer_schrodinger_files.def
```

The build recipe is kept deliberately minimal: it downloads Elmer, copies the
solver source in, and compiles. It does no case setup, which is done in
`cases/`.

## The potential interface

Elmer's `MATC` expressions are used to define the potential. Coordinates are
`tx` (=x), `tx(1)` (=y) and `tx(2)` (=z). Commas separate arguments. Examples:

```
Potential = Variable Coordinate
  Real MATC "1.0*tx*tx"                                  ! 1D harmonic
  Real MATC "0.25*(tx*tx + tx(1)*tx(1))"                 ! 2D harmonic
  Real MATC "0.25*(tx*tx + tx(1)*tx(1) + tx(2)*tx(2))"   ! 3D harmonic
  Real MATC "10.0*(exp(-2*0.5*(tx-3.0))-2*exp(-0.5*(tx-3.0)))"  ! Morse
```

**In two and three dimensions the potential must be separable**,
`V = sum_k V_k(x_k)`. See `docs/MATC_POTENTIAL_SYNTAX.md`.
