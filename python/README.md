# Python reference implementation

An independent implementation of the same discrete problem, written from
scratch, used for code-to-code verification against the Elmer module. The two
codes share no source code, use different mesh containers and different
assembly loops, and solve with different libraries (LAPACK/ARPACK against
Elmer's own solvers).

## Files

| file | purpose |
|---|---|
| `fem_schrodinger.py` | mesh generation, assembly, dense solve |
| `eigsolver_esparso.py` | sparse shift-and-invert Arnoldi solve (ARPACK) |
| `paralelo/` | matrix export to Matrix Market and a PETSc/SLEPc driver |
| `casos/` | one folder per benchmark case, with `run.py` and `command.txt` |

## Running

```bash
python3 fem_schrodinger.py --case casos/1D_harmonico/run.py
```

## Dense against sparse

`fem_schrodinger.py` solves the generalized eigenproblem directly with LAPACK,
which costs `O(n_f^3)` time and `O(n_f^2)` memory. `eigsolver_esparso.py`
factorizes `(A - sigma M)` once and iterates, costing `O(nnz)` storage. On
identical matrices the sparse route is up to 4097 times faster and reaches
problem sizes the dense route cannot.

```python
from fem_schrodinger import malha_2d, montar, contorno
from eigsolver_esparso import resolve_esparso
dims = malha_2d(20.0, 20.0, 200, 200, -10.0, -10.0)
A, M, gmap, X, Vv = montar(dims, lambda X: 0.25*(X[:,0]**2 + X[:,1]**2))
bc = contorno(X, dims)
vals, psi, interno, info = resolve_esparso(A, M, nmodos=4, nos_contorno=bc)
```

## Note on the Dirichlet treatment

Boundary degrees of freedom are **eliminated** (the rows and columns are
removed), not shifted to the diagonal. The diagonal-shift treatment used for
source problems produces spurious eigenvalues here.
