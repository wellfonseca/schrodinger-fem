#!/usr/bin/env python3
"""
Export the FEM matrices (A = K + Ma  and  M) to Matrix Market,
to be read by the parallel PETSc/SLEPc solver.

Usage:  exporta_mtx.py <dimension> <L> <n> [alpha]
      1 / 2 / 3
"""
import sys, os
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
from schrodinger_fem import mesh_1d, mesh_2d, mesh_3d, assemble, boundary


def exporta(dim, L, n, alpha, pasta):
    os.makedirs(pasta, exist_ok=True)
    if dim == 1:
        dims = mesh_1d(L, n, x0=-L/2)
        Vf = lambda X: alpha * X[:, 0]**2
    elif dim == 2:
        dims = mesh_2d(L, L, n, n, -L/2, -L/2)
        Vf = lambda X: alpha * (X[:, 0]**2 + X[:, 1]**2)
    else:
        dims = mesh_3d(L, L, L, n, n, n, -L/2, -L/2, -L/2)
        Vf = lambda X: alpha * (X[:, 0]**2 + X[:, 1]**2 + X[:, 2]**2)

    A, M, gmap, X, V = assemble(dims, Vf)
    bc = boundary(X, dims)
    interior = np.setdiff1d(np.arange(A.shape[0]), bc)

    A = A.tocsr()[interior][:, interior].tocsr()
    M = M.tocsr()[interior][:, interior].tocsr()

    mmwrite(os.path.join(pasta, 'A.mtx'), csr_matrix(A), field='real', symmetry='symmetric')
    mmwrite(os.path.join(pasta, 'M.mtx'), csr_matrix(M), field='real', symmetry='symmetric')

    n_ = A.shape[0]
    print("  dim=%d  L=%g  n=%d  alpha=%g" % (dim, L, n, alpha))
    print("  free DOF   = %d   nnz(A) = %d   nnz(M) = %d" % (n_, A.nnz, M.nnz))
    print("  density    = %.5f %%" % (100.0 * (A.nnz + M.nnz) / (2.0 * n_ * n_)))
    print("  dense memory (3x n^2 x 8) = %.0f MB" % (24.0 * n_ * n_ / 1e6))
    print("  files in %s" % pasta)
    return n_, A.nnz


if __name__ == '__main__':
    dim = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    L = float(sys.argv[2]) if len(sys.argv) > 2 else 30.0
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    alpha = float(sys.argv[4]) if len(sys.argv) > 4 else 0.25
    pasta = sys.argv[5] if len(sys.argv) > 5 else 'matrices'
    exporta(dim, L, n, alpha, pasta)
