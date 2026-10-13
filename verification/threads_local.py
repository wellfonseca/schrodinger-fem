import os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
from fem_schrodinger import malha_2d, montar, contorno
from eigsolver_esparso import resolve_esparso
L,n=20.0,200
dims=malha_2d(L,L,n,n,-L/2,-L/2)
A,M,g,X,V=montar(dims, lambda X: 0.25*(X[:,0]**2+X[:,1]**2))
bc=contorno(X,dims)
print("  GDL livres = %d" % (A.shape[0]-len(bc)))
t=[]
for rep in range(2):
    r=resolve_esparso(A,M,nmodos=4,nos_contorno=bc,verbose=False)[3]
    t.append(r['tempo_s'])
print("  %.3f s" % min(t))
