#!/usr/bin/env python3
"""
Solver FEM para a equacao de Schrodinger estacionaria.

    -lap(psi) + V(x) psi = E psi          (hbar = m = 1)

Discretizacao: elementos lagrangianos em malha estruturada (1D/2D/3D).
Autoproblema generalizado:  A psi = E M psi,   A = K + Ma

  K  = integral( grad(Ni) . grad(Nj) )      rigidez
  M  = integral( Ni Nj )                    massa
  Ma = integral( Ni V Nj )                  massa ponderada pelo potencial

Matrizes locais (elemento linear, lado h):
  K_loc = (1/h) [[1,-1],[-1,1]]
  M_loc = (h/6) [[2,1],[1,2]]

Referencia: Sehra, "Finite Element Analysis of the Schrodinger Equation",
Eqs. 3.20-3.27 (PDF anexo).
"""
import numpy as np


# ----------------------------------------------------------------------
#  malhas estruturadas
# ----------------------------------------------------------------------
def malha_1d(L, n, x0=0.0):
    """n elementos, n+1 nos."""
    x = np.linspace(x0, x0 + L, n + 1)
    return [x]


def malha_2d(Lx, Ly, nx, ny, x0=0.0, y0=0.0):
    x = np.linspace(x0, x0 + Lx, nx + 1)
    y = np.linspace(y0, y0 + Ly, ny + 1)
    return [x, y]


def malha_3d(Lx, Ly, Lz, nx, ny, nz, x0=0.0, y0=0.0, z0=0.0):
    x = np.linspace(x0, x0 + Lx, nx + 1)
    y = np.linspace(y0, y0 + Ly, ny + 1)
    z = np.linspace(z0, z0 + Lz, nz + 1)
    return [x, y, z]


# ----------------------------------------------------------------------
#  assembleia
# ----------------------------------------------------------------------
def _nos(dims):
    """Numero de nos por dimensao."""
    return [len(a) - 1 for a in dims]     # nº de elementos por direcao


def mapa_nos(dims):
    """Indice global de cada no da malha estruturada."""
    ne = _nos(dims)
    if len(dims) == 1:
        return np.arange(ne[0] + 1)
    if len(dims) == 2:
        nx, ny = ne
        return np.arange((nx + 1) * (ny + 1)).reshape(nx + 1, ny + 1, order='F')
    nx, ny, nz = ne
    return np.arange((nx + 1) * (ny + 1) * (nz + 1)).reshape(
        nx + 1, ny + 1, nz + 1, order='F')


def montar(dims, Vfun, elementos=1):
    """Monta A = K + Ma e M.

    dims     : lista de arrays de coordenadas [x] ou [x,y] ou [x,y,z]
    Vfun     : callable(X) -> V, recebendo array (nnos, ndim)
    elementos: 1 = linear
    """
    ndim = len(dims)
    ne = _nos(dims)
    gmap = mapa_nos(dims)
    nnos = gmap.size

    # coordenadas de todos os nos
    if ndim == 1:
        X = np.column_stack([dims[0]])
    elif ndim == 2:
        Y, X = np.meshgrid(dims[1], dims[0], indexing='ij')
        X = np.column_stack([X.ravel(order='F'), Y.ravel(order='F')])
    else:
        Z, Y, X = np.meshgrid(dims[2], dims[1], dims[0], indexing='ij')
        X = np.column_stack([X.ravel(order='F'), Y.ravel(order='F'), Z.ravel(order='F')])

    V = np.asarray(Vfun(X), dtype=float).ravel()

    from scipy.sparse import lil_matrix
    A = lil_matrix((nnos, nnos))
    M = lil_matrix((nnos, nnos))

    # indices locais dos 2^ndim vertices de um elemento
    loc = np.array(np.meshgrid(*[range(2)] * ndim, indexing='ij')).reshape(ndim, -1).T

    it = np.ndindex(*ne)
    for e in it:
        # indices globais dos nos do elemento
        idx = [gmap[tuple(e[d] + loc[k, d] for d in range(ndim))] for k in range(loc.shape[0])]
        idx = np.array(idx)

        # h por direcao
        h = np.array([dims[d][1] - dims[d][0] for d in range(ndim)])

        # matrizes locais por QUADRATURA DE GAUSS (correta em 1D/2D/3D)
        # (o produto tensorial das matrizes 1D NAO gera o elemento bilinear)
        Kl, Ml, Mal = _locais(ndim, h, V[idx])

        for a in range(len(idx)):
            for b in range(len(idx)):
                A[idx[a], idx[b]] += Kl[a, b] + Mal[a, b]
                M[idx[a], idx[b]] += Ml[a, b]

    return A.tocsr(), M.tocsr(), gmap, X, V


def _locais(ndim, h, Vv):
    """Matrizes locais (rigidez, massa, massa ponderada) por quadratura de Gauss.

    Elemento lagrangiano multilinear. Derivadas:
        dN/dx_d = (2/h_d) * dN/dxi_d
    Integracao: 2 pontos de Gauss por direcao (exata para estes polinomios).
    """
    g = np.array([-1.0, 1.0]) / np.sqrt(3.0)
    w = np.array([1.0, 1.0])
    nloc = 2 ** ndim
    # bits de cada no local: formas[k][d] = 0/1
    formas = [[(k >> d) & 1 for d in range(ndim)] for k in range(nloc)]

    def N(k, xi):
        v = 1.0
        for d in range(ndim):
            v *= (1 + xi[d]) / 2 if formas[k][d] else (1 - xi[d]) / 2
        return v

    def dN(k, d, xi):
        v = (0.5 if formas[k][d] else -0.5)
        for e in range(ndim):
            if e == d:
                continue
            v *= (1 + xi[e]) / 2 if formas[k][e] else (1 - xi[e]) / 2
        return v

    Kl = np.zeros((nloc, nloc))
    Ml = np.zeros((nloc, nloc))
    Mal = np.zeros((nloc, nloc))
    detJ = np.prod(h) / (2 ** ndim)          # |J| = prod(h_d/2)

    for gp in np.ndindex(*([2] * ndim)):
        xi = np.array([g[i] for i in gp])
        peso = np.prod([w[i] for i in gp]) * detJ

        Ni = np.array([N(k, xi) for k in range(nloc)])
        dNi = np.array([[dN(k, d, xi) * (2.0 / h[d]) for d in range(ndim)]
                        for k in range(nloc)])

        Kl += peso * (dNi @ dNi.T)
        Ml += peso * np.outer(Ni, Ni)
        Vg = float(np.dot(Ni, Vv))
        Mal += peso * Vg * np.outer(Ni, Ni)

    return Kl, Ml, Mal


# ----------------------------------------------------------------------
#  solucao do autoproblema
# ----------------------------------------------------------------------
def contorno(X, dims, tol=1e-12):
    """Nos de contorno (Dirichlet psi=0), identificados pela POSICAO.

    X    : array (nnos, ndim) de coordenadas
    dims : lista de arrays de coordenadas
    """
    ndim = len(dims)
    mask = np.zeros(len(X), dtype=bool)
    for d in range(ndim):
        lo, hi = dims[d][0], dims[d][-1]
        mask |= np.isclose(X[:, d], lo, atol=tol) | np.isclose(X[:, d], hi, atol=tol)
    return np.where(mask)[0]


def resolver(A, M, nmodos=5, nos_contorno=None, metodo='auto', denso_max=8000):
    """Resolve A psi = E M psi.

    metodo: 'denso' (LAPACK), 'sparse' (ARPACK shift-invert) ou 'auto'.
    """
    A = A.tocsr(); M = M.tocsr()
    if nos_contorno is not None and len(nos_contorno):
        interno = np.setdiff1d(np.arange(A.shape[0]), nos_contorno)
        A = A[interno][:, interno]
        M = M[interno][:, interno]
    else:
        interno = np.arange(A.shape[0])

    n = A.shape[0]
    if metodo == 'auto':
        metodo = 'denso' if n <= denso_max else 'sparse'

    if metodo == 'denso':
        from scipy.linalg import eigh
        E, vec = eigh(A.toarray(), M.toarray())
        ordem = np.argsort(E)
        E = E[ordem][:nmodos]
        vec = vec[:, ordem][:, :nmodos]
    else:
        from scipy.sparse.linalg import eigsh
        E, vec = eigsh(A, k=nmodos, M=M, sigma=0.0, which='LM')
        ordem = np.argsort(E)
        E = E[ordem]; vec = vec[:, ordem]

    # reconstruir com zeros no contorno
    psi = np.zeros((A.shape[0] + (len(nos_contorno) if nos_contorno is not None else 0), len(E)))
    psi[interno, :] = vec
    return E, psi, interno


def normalizar(psi, M):
    """Normaliza cada modo: psi^T M psi = 1."""
    for j in range(psi.shape[1]):
        nrm = np.sqrt(psi[:, j] @ (M @ psi[:, j]))
        if nrm > 0:
            psi[:, j] /= nrm
    return psi
