#!/usr/bin/env python3
"""
FEM solver for the stationary Schrodinger equation.

    -lap(psi) + V(x) psi = E psi          (hbar = m = 1)

Discretisation: Lagrangian elements on a structured mesh (1D/2D/3D).
Generalised eigenproblem:  A psi = E M psi,   A = K + Ma

  K  = integral( grad(Ni) . grad(Nj) )      stiffness
  M  = integral( Ni Nj )                    mass
  Ma = integral( Ni V Nj )                  potential-weighted mass matrix

Local matrices (linear element, side h):
  K_loc = (1/h) [[1,-1],[-1,1]]
  M_loc = (h/6) [[2,1],[1,2]]

Reference: Sehra, "Finite Element Analysis of the Schrodinger Equation",
Eqs. 3.20-3.27 (attached PDF).
"""
import numpy as np


# ----------------------------------------------------------------------
#  structured meshes
# ----------------------------------------------------------------------
def mesh_1d(L, n, x0=0.0):
    """n elements, n+1 nodes."""
    x = np.linspace(x0, x0 + L, n + 1)
    return [x]


def mesh_2d(Lx, Ly, nx, ny, x0=0.0, y0=0.0):
    x = np.linspace(x0, x0 + Lx, nx + 1)
    y = np.linspace(y0, y0 + Ly, ny + 1)
    return [x, y]


def mesh_3d(Lx, Ly, Lz, nx, ny, nz, x0=0.0, y0=0.0, z0=0.0):
    x = np.linspace(x0, x0 + Lx, nx + 1)
    y = np.linspace(y0, y0 + Ly, ny + 1)
    z = np.linspace(z0, z0 + Lz, nz + 1)
    return [x, y, z]


# ----------------------------------------------------------------------
#  assembly
# ----------------------------------------------------------------------
def _nodes(dims):
    """Number of nodes per dimension."""
    return [len(a) - 1 for a in dims]     # number of elements per direction


def node_map(dims):
    """Global index of each node of the structured mesh."""
    ne = _nodes(dims)
    if len(dims) == 1:
        return np.arange(ne[0] + 1)
    if len(dims) == 2:
        nx, ny = ne
        return np.arange((nx + 1) * (ny + 1)).reshape(nx + 1, ny + 1, order='F')
    nx, ny, nz = ne
    return np.arange((nx + 1) * (ny + 1) * (nz + 1)).reshape(
        nx + 1, ny + 1, nz + 1, order='F')


def assemble(dims, Vfun, n_elements=1):
    """Assemble A = K + Ma and M.

    dims     : list of coordinate arrays [x] or [x,y] or [x,y,z]
    Vfun     : callable(X) -> V, receiving an array (nnos, ndim)
    n_elements: 1 = linear
    """
    ndim = len(dims)
    ne = _nodes(dims)
    gmap = node_map(dims)
    nnos = gmap.size

    # coordinates of all nodes
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

    # local indices of the 2^ndim vertices of an element
    loc = np.array(np.meshgrid(*[range(2)] * ndim, indexing='ij')).reshape(ndim, -1).T

    it = np.ndindex(*ne)
    for e in it:
        # global indices of the element nodes
        idx = [gmap[tuple(e[d] + loc[k, d] for d in range(ndim))] for k in range(loc.shape[0])]
        idx = np.array(idx)

        # h per direction
        h = np.array([dims[d][1] - dims[d][0] for d in range(ndim)])

        # local matrices by GAUSS QUADRATURE (correct in 1D/2D/3D)
        # (the tensor product of the 1D matrices does NOT generate the bilinear element)
        Kl, Ml, Mal = _local_matrices(ndim, h, V[idx])

        for a in range(len(idx)):
            for b in range(len(idx)):
                A[idx[a], idx[b]] += Kl[a, b] + Mal[a, b]
                M[idx[a], idx[b]] += Ml[a, b]

    return A.tocsr(), M.tocsr(), gmap, X, V


def _local_matrices(ndim, h, Vv):
    """Local matrices (stiffness, mass, potential-weighted mass) by Gauss quadrature.

    Multilinear Lagrangian element. Derivatives:
        dN/dx_d = (2/h_d) * dN/dxi_d
    Integration: 2 Gauss points per direction (exact for these polynomials).
    """
    g = np.array([-1.0, 1.0]) / np.sqrt(3.0)
    w = np.array([1.0, 1.0])
    nloc = 2 ** ndim
    # bits of each local node: formas[k][d] = 0/1
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
#  solution of the eigenproblem
# ----------------------------------------------------------------------
def boundary(X, dims, tol=1e-12):
    """Boundary nodes (Dirichlet psi=0), identified by POSITION.

    X    : array (nnos, ndim) of coordinates
    dims : list of coordinate arrays
    """
    ndim = len(dims)
    mask = np.zeros(len(X), dtype=bool)
    for d in range(ndim):
        lo, hi = dims[d][0], dims[d][-1]
        mask |= np.isclose(X[:, d], lo, atol=tol) | np.isclose(X[:, d], hi, atol=tol)
    return np.where(mask)[0]


def solve(A, M, nmodes=5, boundary_nodes=None, method='auto', dense_max=8000):
    """Solve A psi = E M psi.

    method: 'dense' (LAPACK), 'sparse' (ARPACK shift-invert) or 'auto'.
    """
    A = A.tocsr(); M = M.tocsr()
    if boundary_nodes is not None and len(boundary_nodes):
        interior = np.setdiff1d(np.arange(A.shape[0]), boundary_nodes)
        A = A[interior][:, interior]
        M = M[interior][:, interior]
    else:
        interior = np.arange(A.shape[0])

    n = A.shape[0]
    if method == 'auto':
        method = 'dense' if n <= dense_max else 'sparse'

    if method == 'dense':
        from scipy.linalg import eigh
        E, vec = eigh(A.toarray(), M.toarray())
        order = np.argsort(E)
        E = E[order][:nmodes]
        vec = vec[:, order][:, :nmodes]
    else:
        from scipy.sparse.linalg import eigsh
        E, vec = eigsh(A, k=nmodes, M=M, sigma=0.0, which='LM')
        order = np.argsort(E)
        E = E[order]; vec = vec[:, order]

    # rebuild with zeros on the boundary
    psi = np.zeros((A.shape[0] + (len(boundary_nodes) if boundary_nodes is not None else 0), len(E)))
    psi[interior, :] = vec
    return E, psi, interior


def normalize(psi, M):
    """Normalise each mode: psi^T M psi = 1."""
    for j in range(psi.shape[1]):
        nrm = np.sqrt(psi[:, j] @ (M @ psi[:, j]))
        if nrm > 0:
            psi[:, j] /= nrm
    return psi
