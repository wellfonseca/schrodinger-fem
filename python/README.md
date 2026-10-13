# Python reference implementation

An independent implementation of the same discrete problem, written from
scratch, used for code-to-code verification against the Elmer module. The two
codes share no source code, use different mesh containers and different
assembly loops, and solve with different libraries (LAPACK/ARPACK against
Elmer's own solvers).

## Files

| file | purpose |
|---|---|
| `schrodinger_fem.py` | mesh generation, assembly, dense solve |
| `sparse_eigsolver.py` | sparse shift-and-invert Arnoldi solve (ARPACK) |
| `paralelo/` | matrix export to Matrix Market and a PETSc/SLEPc driver |
| `cases/` | one folder per benchmark case, with `run.py` and `command.txt` |

## Running

```bash
python3 cases/1D_harmonic/run.py
```

## Dense against sparse

`schrodinger_fem.py` solves the generalized eigenproblem directly with LAPACK,
which costs `O(n_f^3)` time and `O(n_f^2)` memory. `sparse_eigsolver.py`
factorizes `(A - sigma M)` once and iterates, costing `O(nnz)` storage. On
identical matrices the sparse route is up to 4097 times faster and reaches
problem sizes the dense route cannot.

```python
from schrodinger_fem import mesh_2d, assemble, boundary
from sparse_eigsolver import solve_sparse
dims = mesh_2d(20.0, 20.0, 200, 200, -10.0, -10.0)
A, M, gmap, X, Vv = assemble(dims, lambda X: 0.25*(X[:,0]**2 + X[:,1]**2))
bc = boundary(X, dims)
vals, psi, interior, info = solve_sparse(A, M, nmodes=4, boundary_nodes=bc)
```

## Note on the Dirichlet treatment

Boundary degrees of freedom are **eliminated** (the rows and columns are
removed), not shifted to the diagonal. The diagonal-shift treatment used for
source problems produces spurious eigenvalues here.
