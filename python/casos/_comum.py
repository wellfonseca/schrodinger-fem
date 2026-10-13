#!/usr/bin/env python3
"""Utilitarios comuns aos casos do schrodinger_P1_1."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
from fem_schrodinger import (malha_1d, malha_2d, malha_3d, montar,
                             contorno, resolver, normalizar)


def roda(dim, L, n, alpha, modo, ana_fn, nmodos=5):
    """Monta, resolve e imprime a comparacao com o analitico."""
    print("=" * 64)
    print(" Caso: %s" % modo)
    print("=" * 64)
    print("  %dD   L = %.1f   n = %d   alpha = %.4f" % (dim, L, n, alpha))
    print("  convencao: -lap(psi) + V psi = E psi")
    print()

    if dim == 1:
        dims = malha_1d(L, n, x0=-L/2)
    elif dim == 2:
        dims = malha_2d(L, L, n, n, -L/2, -L/2)
    else:
        dims = malha_3d(L, L, L, n, n, n, -L/2, -L/2, -L/2)

    if alpha == 0.0:
        Vf = lambda X: np.zeros(len(X))
    elif dim == 1:
        Vf = lambda X: alpha * X[:, 0]**2
    elif dim == 2:
        Vf = lambda X: alpha * (X[:, 0]**2 + X[:, 1]**2)
    else:
        Vf = lambda X: alpha * (X[:, 0]**2 + X[:, 1]**2 + X[:, 2]**2)

    A, M, gmap, X, V = montar(dims, Vf)
    bc = contorno(X, dims)
    E, psi, interno = resolver(A, M, nmodos=nmodos, nos_contorno=bc)
    psi = normalizar(psi, M)

    print("  GDL = %d   (contorno %d, livres %d)" % (A.shape[0], len(bc), len(interno)))
    print("  V_min = %.6f   V_max = %.6f" % (V.min(), V.max()))
    print()
    print("  %-4s %-18s %-18s %-12s" % ("n", "E calculado", "E analitico", "erro rel."))
    for k, ek in enumerate(E):
        a = ana_fn(k)
        if a is None:
            print("  %-4d %-18.10f %-18s %-12s" % (k + 1, ek, "(n/d)", "-"))
        else:
            print("  %-4d %-18.10f %-18.10f %-12.3e"
                  % (k + 1, ek, a, abs(ek - a) / abs(a)))

    nrm = np.array([psi[:, j] @ (M @ psi[:, j]) for j in range(psi.shape[1])])
    print()
    print("  norma psi^T M psi: min = %.10f  max = %.10f  (deve ser 1)"
          % (nrm.min(), nrm.max()))
    return E, A, M, psi
