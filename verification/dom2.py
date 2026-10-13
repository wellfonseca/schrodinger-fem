import numpy as np
def E1_2d(L,h=0.25):
    n=int(round(L/h)); N=n+1; x=np.linspace(-L/2,L/2,N)
    g2=np.array([-1/np.sqrt(3),1/np.sqrt(3)]); w2=np.array([1.0,1.0])
    K=np.zeros((N,N));M=np.zeros((N,N));Vi=np.zeros((N,N))
    for i in range(n):
        a,b=x[i],x[i+1]; hh=b-a; J=hh/2; Vn=np.array([a*a,b*b])
        for xi,wi in zip(g2,w2):
            Nv=np.array([(1-xi)/2,(1+xi)/2]); dN=np.array([-0.5,0.5])/J
            M[i:i+2,i:i+2]+=wi*J*np.outer(Nv,Nv); K[i:i+2,i:i+2]+=wi*J*np.outer(dN,dN)
            Vi[i:i+2,i:i+2]+=wi*J*float(np.dot(Nv,Vn))*np.outer(Nv,Nv)
    A=np.kron(K,M)+np.kron(M,K)+0.25*(np.kron(Vi,M)+np.kron(M,Vi)); Mm=np.kron(M,M)
    livre=np.ones((N,N),bool); livre[0]=livre[-1]=False; livre[:,0]=livre[:,-1]=False
    f=np.where(livre.ravel())[0]
    Lc=np.linalg.cholesky(Mm[np.ix_(f,f)])
    C=np.linalg.solve(Lc,A[np.ix_(f,f)]); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
    return np.linalg.eigvalsh(C)[0], len(f)
print("═══ TAMANHO DE DOMINIO a MALHA FIXA (h=0.25) — 2D harmonico, exato 1.0 ═══")
print("     L      n_f      E1                E1 - 1.0")
for L in [12,16,20,24,28]:
    e,nf=E1_2d(L,0.25)
    print("  %5.1f  %7d  %.12f   %+.3e" % (L,nf,e,e-1.0))
