# schrodinger-fem v1.1.0

Adds the InGaN/GaN benchmark requested in review, and position-dependent
effective mass.

## Position-dependent effective mass

The solver now assembles the conservative form
`-(hbar^2/2) div[(1/m*) grad psi] + V psi = E psi`, which is the correct
Hamiltonian for a heterostructure. The stiffness carries a factor `1/m*`
evaluated at the Gauss point, read from a new optional Body Force keyword:

```
Effective Mass = Variable Coordinate 1
  Real MATC "0.190 - 0.012*( (tx>=3.53)*(tx<=6.53) + (tx>=9.53)*(tx<=12.53) )"
```

The default is `m* = 1`, so every existing case is unchanged bit for bit:
`1D_harmonic` still gives `E1 = 1.000572892257`.

Validated against the independent Python implementation on the InGaN/GaN
double well of En-nadir et al. (2022), where the mass goes from 0.178 m0 in
the wells to 0.190 m0 in the barriers:

| mass model | Elmer (meV) | Python (meV) |
|---|---|---|
| constant 0.190 m0 | 112.838 / 113.786 / 414.684 | 112.83 / 113.78 / 414.66 |
| variable | 116.192 / 117.216 / 427.171 | 116.19 / 117.20 / 427.06 |

A 6.3 % mass contrast moves the levels by 3.0 %, two orders of magnitude
above the discretization error of the case.

## New benchmark case

`cases/1D_ingan_gan_dqw/` — symmetric InGaN/GaN double quantum well after
En-nadir et al., *J. Comput. Electron.* **21** (2022) 111-118,
[doi:10.1007/s10825-021-01830-4](https://doi.org/10.1007/s10825-021-01830-4).

Two 3 nm InGaN wells separated by a 3 nm GaN barrier, `V0 = 0.616 eV` from
`Eg(InN) = 0.78 eV`, `Eg(GaN) = 3.42 eV`, bowing `b = 3.8 eV`, `x_In = 0.15`
and a 70/30 band offset.

The well interfaces are deliberately placed off the mesh nodes. When an
interface coincides with a node, the observed convergence order drops from
two to one; the two sequences converge to the same limit, so no spurious
value is produced, but the rate halves. This is documented in the case's
`command.txt` and in the manuscript.

## What is still not covered

- The donor impurity and the dielectric mismatch of the reference are not
  represented; only the structure, the spectrum and the field dependence are.
- The flux continuity `(1/m*) dpsi/dn` across an interface is satisfied only
  in the weak sense, because `V_h` and `m*_h` are nodal interpolants of
  piecewise-constant fields.
- Self-consistent Schrodinger-Poisson is outside the present scope.

## Installing

The module is a single Fortran source. Compile against an Elmer 9.0 build:

```
gfortran -c SchrodingerSolver.F90 -o SchrodingerSolver.o \
  -DCONTIG= -DHAVE_EXECUTECOMMANDLINE -DUSE_ISO_C_BINDINGS \
  -O2 -fPIC -I<elmer-build>/fmodules
gfortran -shared SchrodingerSolver.o -o SchrodingerSolver.so \
  -L<elmer-install>/lib/elmersolver -lelmersolver
```

and put `SchrodingerSolver.so` where Elmer looks for user functions.

## Verification

`verification/smoke_test.sh` runs every Python case end to end.
`verification/timings/README.md` records the raw timing data behind the
performance table, from two measurement sessions.
