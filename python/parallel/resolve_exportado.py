#!/usr/bin/env python3
"""
Solve A x = E M x from the matrices exported by the SchrodingerSolver.

Flow (the solution for parallelism, since Elmer 9.0 does not partition a mesh
on this platform and the solver assembles dense local matrices):

    Elmer (assembles A and M, exports)  ->  this script (solves SPARSE)

Advantages:
  - uses ARPACK shift-invert (sparse): memory O(nnz), not O(n^2)
  - measured speed-up: up to 4097x vs dense
  - independent of Elmer's ARPACK and of the dense DSYGV

Usage:
    python3 resolve_exportado.py <prefix> [nmodes]
    python3 resolve_exportado.py schrod 5
"""
import sys, os, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def le_matriz(caminho):
    """Reads the ExportaMatriz format: '% ...' , 'n n nnz', then 'i j v' (0-based)."""
    with open(caminho) as f:
        corpo = [l for l in f if not l.startswith('%')]
    n, m, nnz_declarado = map(int, corpo[0].split())
    linhas, cols, vals = [], [], []
    for l in corpo[1:]:
        p = l.split()
        if len(p) >= 3:
            linhas.append(int(p[0])); cols.append(int(p[1])); vals.append(float(p[2]))
    M = coo_matrix((vals, (linhas, cols)), shape=(n, m)).tocsr()
    return M, n, nnz_declarado


def le_fixos(caminho):
    with open(caminho) as f:
        linhas = [l.strip() for l in f if l.strip() != '']
    nfix = int(linhas[0])
    return np.array([int(x) for x in linhas[1:1+nfix]], dtype=int)


def main():
    prefix = sys.argv[1] if len(sys.argv) > 1 else 'schrod'
    nmodes = int(sys.argv[2]) if len(sys.argv) > 2 else 5

    print("=" * 64)
    print(" Solving the system exported by the SchrodingerSolver")
    print("=" * 64)
    print("  prefix: %s" % prefix)

    t0 = time.time()
    A, n, nnzA = le_matriz('%s_A.dat' % prefix)
    M, n2, nnzM = le_matriz('%s_M.dat' % prefix)
    assert n == n2, "dimensions of A and M are incompatible"

    print("  A: %d x %d, nnz = %d" % (A.shape[0], A.shape[1], A.nnz))
    print("  M: %d x %d, nnz = %d" % (M.shape[0], M.shape[1], M.nnz))
    dens_densa = 24.0 * n * n / 1e6
    dens_esparsa = 12.0 * (A.nnz + M.nnz) / 1e6
    print("  memory: dense %.1f MB  |  sparse ~%.1f MB  (%.0fx smaller)"
          % (dens_densa, dens_esparsa, dens_densa / max(dens_esparsa, 1e-9)))

    # eliminate the Dirichlet DOFs
    if os.path.exists('%s_fixed.dat' % prefix):
        fixos = le_fixos('%s_fixed.dat' % prefix)
        interior = np.setdiff1d(np.arange(n), fixos)
        print("  fixed DOFs (Dirichlet): %d  ->  free: %d" % (len(fixos), len(interior)))
        A = A[interior][:, interior].tocsr()
        M = M[interior][:, interior].tocsr()

    print()
    print("  solving %d eigenvalues (ARPACK shift-invert)..." % nmodes)
    t1 = time.time()
    vals, vecs = eigsh(A, k=nmodes, M=M, sigma=0.0, which='LM', tol=0.0)
    dt = time.time() - t1
    order = np.argsort(vals)
    vals = vals[order]

    print("  time: %.2f s   (total %.2f s)" % (dt, time.time() - t0))
    print()
    print("  --- EIGENVALUES ---")
    for i, v in enumerate(vals):
        print("    E(%d) = %+.12e" % (i + 1, v))


if __name__ == '__main__':
    main()
