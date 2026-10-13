import numpy as np
L,n,alpha=20.0,20,0.25
g2=np.array([-1/np.sqrt(3),1/np.sqrt(3)]); w2=np.array([1.0,1.0])
N=n+1; x=np.linspace(-L/2,L/2,N)
K=np.zeros((N,N));M=np.zeros((N,N));Vi=np.zeros((N,N))
for i in range(n):
    a,b=x[i],x[i+1]; h=b-a; J=h/2; Vn=np.array([a*a,b*b])
    for xi,wi in zip(g2,w2):
        Nv=np.array([(1-xi)/2,(1+xi)/2]); dN=np.array([-0.5,0.5])/J
        M[i:i+2,i:i+2]+=wi*J*np.outer(Nv,Nv); K[i:i+2,i:i+2]+=wi*J*np.outer(dN,dN)
        Vi[i:i+2,i:i+2]+=wi*J*float(np.dot(Nv,Vn))*np.outer(Nv,Nv)
A=np.kron(K,M)+np.kron(M,K)+alpha*(np.kron(Vi,M)+np.kron(M,Vi))
Mm=np.kron(M,M)
livre=np.ones((N,N),bool); livre[0]=livre[-1]=False; livre[:,0]=livre[:,-1]=False
f=np.where(livre.ravel())[0]
Af=A[np.ix_(f,f)]; Mf=Mm[np.ix_(f,f)]
print("2D harmonic, L=20, n=20 (free DOF = %d)" % len(f))
Lc=np.linalg.cholesky(Mf); C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
v=np.linalg.eigvalsh(C)[:4]
art=[1.113644656498, 2.172208383659, 2.172208383659, 3.230772110820]
print()
print("  Table 6 of the paper - 'in-framework (dense DSYGV)'")
for i,(a,b) in enumerate(zip(art,v),1):
    print("    E%d  paper %.12f   mine %.12f   dif %+.2e" % (i,a,b,a-b))
print()
print("  degeneracy E2-E3: mine %.2e" % (v[2]-v[1]))
