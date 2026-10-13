# Verification scripts

Scripts that reproduce the numerical tables of the manuscript. They use the
Python reference implementation in `../python/` and require `numpy`; some also
require `scipy`.

| script | reproduces |
|---|---|
| `confere.py` | analytic reference spectra against computed values |
| `conv1d.py` | 1D convergence table (eigenvalues) |
| `conv2d_int.py` | 2D convergence, interpolated-potential convention |
| `conv3d_int.py` | 3D convergence |
| `autof3.py` | eigenfunction error in $L^2$ and $H^1$ |
| `dom2.py` | domain-size study at fixed mesh size |
| `quad.py` | quadrature convention: interpolated against exact potential |
| `orto.py` | mass-orthonormality and degeneracies |
| `virial4.py` | virial relation $\langle T\rangle = \langle V\rangle$ |
| `resid2.py` | normalized residuals of the computed eigenpairs |
| `tab6.py` | exported against in-framework eigenvalues |
| `threads_local.py` | thread scaling of the sparse solve |

## Note

These scripts were used to audit the manuscript and several of them exist
because an early version of a table did **not** survive the check. In
particular:

* `dom2.py` separates truncation from discretization by varying the box at
  fixed mesh size; the two effects had been conflated.
* `virial4.py` recomputes the virial ratio from the assembled matrices. An
  earlier standalone script did not normalize the eigenvector and shifted the
  gradient by one sample, which produced a wrong ratio.
* `resid2.py` uses the correct Cholesky back-transformation $u = L^{-\top}y$
  for the generalized eigenproblem; using the eigenvector of the transformed
  problem directly gives a wrong, and plausible-looking, answer.

The scripts print the values; comparing them with the manuscript tables is the
test.
