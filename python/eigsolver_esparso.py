#!/usr/bin/env python3
"""
Autossolver ESPARSO para o problema generalizado  A x = E M x.

Motivacao (item 3): o solver denso (DSYGV / scipy.linalg.eigh) tem
  - memoria  O(n^2)
  - custo    O(n^3)
e estoura acima de ~8000 GDL: medido nesta sessao, n=14641 exigiu 5.1 GB.

Este modulo usa ARPACK em modo SHIFT-INVERT:
  - fatora (A - sigma M) uma unica vez  (esparso, UMFPACK/SuperLU)
  - itera Arnoldi para os k menores autovalores
  - memoria  ~ O(nnz) + O(n.k)          -> MUITO menor que O(n^2)
  - custo    ~ O(k . fatoracao)

O ARPACK (scipy.sparse.linalg.eigsh) foi verificado nesta sessao como correto
para este problema -- ao contrario do wrapper interno do Elmer (EigenSolve.F90),
que tem RESID nao inicializado e devolve info=-9.

Interface compativel com fem_schrodinger.resolver().
"""
import numpy as np


def conta_nnz(A, M=None):
    nnz = A.nnz
    if M is not None:
        nnz += M.nnz
    return nnz


def estima_memoria_mb(n, nnz, k):
    """Memoria aproximada em MB.

    Denso (DSYGV): 3 matrizes n x n em double = 24 n^2 bytes
    Esparso: valores + indices (~12 bytes por nnz) + vetores de Arnoldi (n.k)
    """
    denso = 24.0 * n * n / 1e6
    esparso = (12.0 * nnz + 8.0 * n * (k + 3) * 4) / 1e6
    return denso, esparso


def resolve_esparso(A, M, nmodos=5, nos_contorno=None, sigma=0.0,
                    tol=0.0, maxiter=None, verbose=True):
    """Resolve A x = E M x pelos k MENORES autovalores, via ARPACK shift-invert.

    A, M         : scipy.sparse (CSR)
    nos_contorno : indices a eliminar (Dirichlet psi=0)
    sigma        : shift (0 = menor magnitude)
    """
    from scipy.sparse.linalg import eigsh
    import time, resource

    A = A.tocsr(); M = M.tocsr()
    n_full = A.shape[0]

    if nos_contorno is not None and len(nos_contorno):
        interno = np.setdiff1d(np.arange(n_full), nos_contorno)
        A = A[interno][:, interno].tocsr()
        M = M[interno][:, interno].tocsr()
    else:
        interno = np.arange(n_full)

    n = A.shape[0]
    nnz = conta_nnz(A, M)
    denso_mb, esparso_mb = estima_memoria_mb(n, nnz, nmodos)

    if verbose:
        print("    n=%d  nnz=%d (%.4f%% denso)" % (n, nnz, 100.0*nnz/(n*n)))
        print("    memoria estimada: denso %.1f MB  |  esparso %.1f MB  (%.0fx menor)"
              % (denso_mb, esparso_mb, denso_mb/max(esparso_mb, 1e-9)))

    r0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    t0 = time.time()
    vals, vecs = eigsh(A, k=nmodos, M=M, sigma=sigma, which='LM',
                       tol=tol, maxiter=maxiter)
    dt = time.time() - t0
    r1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

    ordem = np.argsort(vals)
    vals = vals[ordem]; vecs = vecs[:, ordem]

    if verbose:
        print("    tempo ARPACK: %.2f s   |  pico RSS: %.1f MB (delta %.1f MB)"
              % (dt, r1, r1 - r0))

    psi = np.zeros((n_full, len(vals)))
    psi[interno, :] = vecs
    return vals, psi, interno, {"tempo_s": dt, "nnz": nnz,
                                "mem_denso_mb": denso_mb,
                                "mem_esparso_mb": esparso_mb,
                                "rss_mb": r1}


def resolve_denso(A, M, nmodos=5, nos_contorno=None, verbose=True):
    """Referencia: solver denso (mesma matematica, O(n^2)/O(n^3))."""
    from scipy.linalg import eigh
    import time, resource

    A = A.tocsr(); M = M.tocsr()
    n_full = A.shape[0]
    if nos_contorno is not None and len(nos_contorno):
        interno = np.setdiff1d(np.arange(n_full), nos_contorno)
        A = A[interno][:, interno]; M = M[interno][:, interno]
    else:
        interno = np.arange(n_full)

    n = A.shape[0]
    Ad = A.toarray(); Md = M.toarray()
    r0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    t0 = time.time()
    vals, vecs = eigh(Ad, Md)
    dt = time.time() - t0
    r1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

    vals = vals[:nmodos]; vecs = vecs[:, :nmodos]
    if verbose:
        print("    tempo denso: %.2f s   |  pico RSS: %.1f MB (delta %.1f MB)"
              % (dt, r1, r1 - r0))
    psi = np.zeros((n_full, len(vals)))
    psi[interno, :] = vecs
    return vals, psi, interno, {"tempo_s": dt, "rss_mb": r1}
