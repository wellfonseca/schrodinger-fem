import numpy as np
L,n,alpha=20.0,40,0.25
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
Lc=np.linalg.cholesky(Mf); C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
val,y=np.linalg.eigh(C)
U=np.linalg.solve(Lc.T,y)              # back-transformacao correta
G=U.T@(Mf@U)                           # deve ser a identidade
err=np.abs(G-np.eye(G.shape[0])).max()
print("2D harmonico (L=20, n=40): ortonormalidade M")
print("  max |psi_i^T M psi_j - delta_ij| = %.3e" % err)
print("  4 primeiros autovalores: "+" ".join("%.10f"%v for v in val[:4]))
print()
print("  artigo: 'the deviation from the identity is below 1e-15 in the 2D case'")
