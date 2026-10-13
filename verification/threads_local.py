import os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
from schrodinger_fem import mesh_2d, assemble, boundary
from sparse_eigsolver import solve_sparse
L,n=20.0,200
dims=mesh_2d(L,L,n,n,-L/2,-L/2)
A,M,g,X,V=assemble(dims, lambda X: 0.25*(X[:,0]**2+X[:,1]**2))
bc=boundary(X,dims)
print("  free DOF = %d" % (A.shape[0]-len(bc)))
t=[]
for rep in range(2):
    r=solve_sparse(A,M,nmodes=4,boundary_nodes=bc,verbose=False)[3]
    t.append(r['tempo_s'])
print("  %.3f s" % min(t))
