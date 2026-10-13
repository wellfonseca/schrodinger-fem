#!/usr/bin/env python3
"""
Generate a structured mesh in the Elmer format (1D / 2D / 3D).

Produz mesh.header, mesh.nodes, mesh.elements, mesh.boundary.

Codigos de elemento (elements.def do Elmer):
    202  = 2 nodes, linha          (1D)
    404  = 4 nodes, quadrilatero   (2D e superficie 3D)
    808  = 8 nodes, hexaedro       (3D)

Node ordering (checked against meshes produced by ElmerGrid):
    2D quad : (i,j) (i+1,j) (i+1,j+1) (i,j+1)
    3D hex  : (i,j,k) (i+1,j,k) (i+1,j+1,k) (i,j+1,k)
              (i,j,k+1) (i+1,j,k+1) (i+1,j+1,k+1) (i,j+1,k+1)

In mesh.nodes, column 2 is the corner node (-1 by default) and
columns 3..5 are x, y, z.
"""
import os


def _escreve(caminho, linhas):
    with open(caminho, 'w') as f:
        f.write('\n'.join(linhas) + '\n')


def generate_1d(folder, L, n, x0=0.0):
    h = L / n
    nodes = [(x0 + i * h, 0.0, 0.0) for i in range(n + 1)]
    elems = [(i + 1, i + 2) for i in range(n)]
    bnds = [(1, 1, 1, 0, 101, 1), (2, 2, n, 0, 101, n + 1)]
    _write(folder, nodes, elems, bnds, 202, 101)


def generate_2d(folder, Lx, Ly, nx, ny, x0=0.0, y0=0.0):
    hx, hy = Lx / nx, Ly / ny
    nodes = []
    idx = {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            idx[(i, j)] = len(nodes) + 1
            nodes.append((x0 + i * hx, y0 + j * hy, 0.0))
    elems = []
    for j in range(ny):
        for i in range(nx):
            elems.append((idx[(i, j)], idx[(i + 1, j)],
                         idx[(i + 1, j + 1)], idx[(i, j + 1)]))
    # contorno: 4 lados (quad 404 de 2 nodes)
    bnds, t = [], 1
    for i in range(nx):                    # y = y0
        bnds.append((1, t, i + 1, 0, 202, idx[(i, 0)], idx[(i + 1, 0)])); t += 1
    for j in range(ny):                    # x = x0+Lx
        bnds.append((2, t, j + 1, 0, 202, idx[(nx, j)], idx[(nx, j + 1)])); t += 1
    for i in range(nx):                    # y = y0+Ly
        bnds.append((3, t, i + 1, 0, 202, idx[(i + 1, ny)], idx[(i, ny)])); t += 1
    for j in range(ny):                    # x = x0
        bnds.append((4, t, j + 1, 0, 202, idx[(0, j + 1)], idx[(0, j)])); t += 1
    _write(folder, nodes, elems, bnds, 404, 202)


def generate_3d(folder, Lx, Ly, Lz, nx, ny, nz, x0=0.0, y0=0.0, z0=0.0):
    hx, hy, hz = Lx / nx, Ly / ny, Lz / nz
    nodes, idx = [], {}
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                idx[(i, j, k)] = len(nodes) + 1
                nodes.append((x0 + i * hx, y0 + j * hy, z0 + k * hz))
    elems = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                elems.append((idx[(i, j, k)], idx[(i + 1, j, k)],
                             idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)],
                             idx[(i, j, k + 1)], idx[(i + 1, j, k + 1)],
                             idx[(i + 1, j + 1, k + 1)], idx[(i, j + 1, k + 1)]))
    bnds, t = [], 1
    # z = z0 e z = z0+Lz  (faces quad)
    for j in range(ny):
        for i in range(nx):
            bnds.append((1, t, j * nx + i + 1, 0, 404,
                        idx[(i, j, 0)], idx[(i + 1, j, 0)],
                        idx[(i + 1, j + 1, 0)], idx[(i, j + 1, 0)])); t += 1
    for j in range(ny):
        for i in range(nx):
            bnds.append((2, t, j * nx + i + 1, 0, 404,
                        idx[(i, j, nz)], idx[(i, j + 1, nz)],
                        idx[(i + 1, j + 1, nz)], idx[(i + 1, j, nz)])); t += 1
    # y = y0 e y = y0+Ly
    for k in range(nz):
        for i in range(nx):
            bnds.append((3, t, k * nx + i + 1, 0, 404,
                        idx[(i, 0, k)], idx[(i, 0, k + 1)],
                        idx[(i + 1, 0, k + 1)], idx[(i + 1, 0, k)])); t += 1
    for k in range(nz):
        for i in range(nx):
            bnds.append((4, t, k * nx + i + 1, 0, 404,
                        idx[(i, ny, k)], idx[(i + 1, ny, k)],
                        idx[(i + 1, ny, k + 1)], idx[(i, ny, k + 1)])); t += 1
    # x = x0 e x = x0+Lx
    for k in range(nz):
        for j in range(ny):
            bnds.append((5, t, k * ny + j + 1, 0, 404,
                        idx[(0, j, k)], idx[(0, j + 1, k)],
                        idx[(0, j + 1, k + 1)], idx[(0, j, k + 1)])); t += 1
    for k in range(nz):
        for j in range(ny):
            bnds.append((6, t, k * ny + j + 1, 0, 404,
                        idx[(nx, j, k)], idx[(nx, j, k + 1)],
                        idx[(nx, j + 1, k + 1)], idx[(nx, j + 1, k)])); t += 1
    _write(folder, nodes, elems, bnds, 808, 404)


def bnd_bulk_codes(bnds):
    """Codigos de tipo de contorno presentes."""
    return sorted(set(b[4] for b in bnds))


def _write(folder, nodes, elems, bnds, bulk_code, bnd_code):
    os.makedirs(folder, exist_ok=True)
    nb = len(bnds)
    # segunda secao do header: codigo do elemento + No de NOS por elemento
    nnos_bnd = len(bnds[0]) - 5 if bnds else 0
    L = ['%d %d %d' % (len(nodes), len(elems), nb), '%d' % len(bnd_bulk_codes(bnds))]
    for code in sorted(set(b[4] for b in bnds)):
        L.append('%d %d' % (code, nnos_bnd))
    _escreve(os.path.join(folder, 'mesh.header'), L)

    L = []
    for i, (x, y, z) in enumerate(nodes, 1):
        L.append('%d -1 %.10g %.10g %.10g' % (i, x, y, z))
    _escreve(os.path.join(folder, 'mesh.nodes'), L)

    L = []
    for n, e in enumerate(elems, 1):
        L.append('%d 1 %d %s' % (n, bulk_code, ' '.join(str(v) for v in e)))
    _escreve(os.path.join(folder, 'mesh.elements'), L)

    L = []
    for b in bnds:
        L.append('%d %d %d %d %d %s' % (b[0], b[1], b[2], b[3], b[4],
                                        ' '.join(str(v) for v in b[5:])))
    _escreve(os.path.join(folder, 'mesh.boundary'), L)
    print('  %s : %d nodes, %d elementos, %d faces de contorno'
          % (folder, len(nodes), len(elems), nb))


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'teste':
        generate_1d('/tmp/t1', 40.0, 250, -20.0)
        generate_2d('/tmp/t2', 10.0, 10.0, 40, 40)
        generate_3d('/tmp/t3', 6.0, 6.0, 6.0, 10, 10, 10)
