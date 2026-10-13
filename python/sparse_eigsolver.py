#!/usr/bin/env python3
"""
SPARSE eigensolver for the generalised problem  A x = E M x.

Motivation (item 3): the dense solver (DSYGV / scipy.linalg.eigh) has
  - memory  O(n^2)
  - cost    O(n^3)
and overflows above ~8000 DOF: measured in this session, n=14641 required 5.1 GB.

This module uses ARPACK in SHIFT-INVERT mode:
  - factorises (A - sigma M) only once  (sparse, UMFPACK/SuperLU)
  - runs Arnoldi iteration for the k smallest eigenvalues
  - memory  ~ O(nnz) + O(n.k)          -> MUCH smaller than O(n^2)
  - cost    ~ O(k . factorization)

ARPACK (scipy.sparse.linalg.eigsh) was verified in this session as correct
for this problem -- unlike the Elmer interior wrapper (EigenSolve.F90),
which has an uninitialised RESID and returns info=-9.

Interface compatible with schrodinger_fem.solve().
"""
import numpy as np


def count_nnz(A, M=None):
    nnz = A.nnz
    if M is not None:
        nnz += M.nnz
    return nnz


def estimate_memory_mb(n, nnz, k):
    """Approximate memory in MB.

    Dense (DSYGV): 3 n x n matrices in double = 24 n^2 bytes
    Sparse: values + indices (~12 bytes per nnz) + Arnoldi vectors (n.k)
    """
    denso = 24.0 * n * n / 1e6
    esparso = (12.0 * nnz + 8.0 * n * (k + 3) * 4) / 1e6
    return denso, esparso


def solve_sparse(A, M, nmodes=5, boundary_nodes=None, sigma=0.0,
                    tol=0.0, maxiter=None, verbose=True):
    """Solve A x = E M x for the k SMALLEST eigenvalues, via ARPACK shift-invert.

    A, M         : scipy.sparse (CSR)
    boundary_nodes : indices to eliminate (Dirichlet psi=0)
    sigma        : shift (0 = smallest magnitude)
    """
    from scipy.sparse.linalg import eigsh
    import time, resource

    A = A.tocsr(); M = M.tocsr()
    n_full = A.shape[0]

    if boundary_nodes is not None and len(boundary_nodes):
        interior = np.setdiff1d(np.arange(n_full), boundary_nodes)
        A = A[interior][:, interior].tocsr()
        M = M[interior][:, interior].tocsr()
    else:
        interior = np.arange(n_full)

    n = A.shape[0]
    nnz = count_nnz(A, M)
    denso_mb, esparso_mb = estimate_memory_mb(n, nnz, nmodes)

    if verbose:
        print("    n=%d  nnz=%d (%.4f%% dense)" % (n, nnz, 100.0*nnz/(n*n)))
        print("    estimated memory: dense %.1f MB  |  sparse %.1f MB  (%.0fx smaller)"
              % (denso_mb, esparso_mb, denso_mb/max(esparso_mb, 1e-9)))

    r0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    t0 = time.time()
    vals, vecs = eigsh(A, k=nmodes, M=M, sigma=sigma, which='LM',
                       tol=tol, maxiter=maxiter)
    dt = time.time() - t0
    r1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

    order = np.argsort(vals)
    vals = vals[order]; vecs = vecs[:, order]

    if verbose:
        print("    ARPACK time: %.2f s   |  peak RSS: %.1f MB (delta %.1f MB)"
              % (dt, r1, r1 - r0))

    psi = np.zeros((n_full, len(vals)))
    psi[interior, :] = vecs
    return vals, psi, interior, {"tempo_s": dt, "nnz": nnz,
                                "mem_denso_mb": denso_mb,
                                "mem_esparso_mb": esparso_mb,
                                "rss_mb": r1}


def solve_dense(A, M, nmodes=5, boundary_nodes=None, verbose=True):
    """Reference: dense solver (same mathematics, O(n^2)/O(n^3))."""
    from scipy.linalg import eigh
    import time, resource

    A = A.tocsr(); M = M.tocsr()
    n_full = A.shape[0]
    if boundary_nodes is not None and len(boundary_nodes):
        interior = np.setdiff1d(np.arange(n_full), boundary_nodes)
        A = A[interior][:, interior]; M = M[interior][:, interior]
    else:
        interior = np.arange(n_full)

    n = A.shape[0]
    Ad = A.toarray(); Md = M.toarray()
    r0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    t0 = time.time()
    vals, vecs = eigh(Ad, Md)
    dt = time.time() - t0
    r1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

    vals = vals[:nmodes]; vecs = vecs[:, :nmodes]
    if verbose:
        print("    dense time: %.2f s   |  peak RSS: %.1f MB (delta %.1f MB)"
              % (dt, r1, r1 - r0))
    psi = np.zeros((n_full, len(vals)))
    psi[interior, :] = vecs
    return vals, psi, interior, {"tempo_s": dt, "rss_mb": r1}
