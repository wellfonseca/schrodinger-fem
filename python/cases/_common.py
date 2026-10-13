#!/usr/bin/env python3
"""Shared helpers for the Python reference cases."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from schrodinger_fem import (mesh_1d, mesh_2d, mesh_3d, assemble,
                             boundary, solve, normalize)


def run_case(dim, L, n, alpha, label, ana_fn, nmodes=5):
    """Assemble, solve, and print the comparison against the analytic values."""
    print("=" * 64)
    print(" Case: %s" % label)
    print("=" * 64)
    print("  %dD   L = %.1f   n = %d   alpha = %.4f" % (dim, L, n, alpha))
    print("  convention: -lap(psi) + V psi = E psi")
    print()

    if dim == 1:
        dims = mesh_1d(L, n, x0=-L/2)
    elif dim == 2:
        dims = mesh_2d(L, L, n, n, -L/2, -L/2)
    else:
        dims = mesh_3d(L, L, L, n, n, n, -L/2, -L/2, -L/2)

    if alpha == 0.0:
        Vf = lambda X: np.zeros(len(X))
    elif dim == 1:
        Vf = lambda X: alpha * X[:, 0]**2
    elif dim == 2:
        Vf = lambda X: alpha * (X[:, 0]**2 + X[:, 1]**2)
    else:
        Vf = lambda X: alpha * (X[:, 0]**2 + X[:, 1]**2 + X[:, 2]**2)

    A, M, gmap, X, V = assemble(dims, Vf)
    bc = boundary(X, dims)
    E, psi, interior = solve(A, M, nmodes=nmodes, boundary_nodes=bc)
    psi = normalize(psi, M)

    print("  DOFs = %d   (boundary %d, free %d)" % (A.shape[0], len(bc), len(interior)))
    print("  V_min = %.6f   V_max = %.6f" % (V.min(), V.max()))
    print()
    print("  %-4s %-18s %-18s %-12s" % ("n", "computed E", "analytic E", "rel. error"))
    for k, ek in enumerate(E):
        a = ana_fn(k)
        if a is None:
            print("  %-4d %-18.10f %-18s %-12s" % (k + 1, ek, "(n/a)", "-"))
        else:
            print("  %-4d %-18.10f %-18.10f %-12.3e"
                  % (k + 1, ek, a, abs(ek - a) / abs(a)))

    nrm = np.array([psi[:, j] @ (M @ psi[:, j]) for j in range(psi.shape[1])])
    print()
    print("  norm psi^T M psi: min = %.10f  max = %.10f  (should be 1)"
          % (nrm.min(), nrm.max()))
    return E, A, M, psi
