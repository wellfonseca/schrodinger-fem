import numpy as np
ALPHA=0.25
def mats1d(a,b):
    h=b-a; J=h/2
    g2=np.array([-1/np.sqrt(3),1/np.sqrt(3)]); w2=np.array([1.0,1.0])
    K=np.zeros((2,2));M=np.zeros((2,2));Vi=np.zeros((2,2))
    Vn=np.array([a*a,b*b])
    for xi,wi in zip(g2,w2):
        Nv=np.array([(1-xi)/2,(1+xi)/2]); dN=np.array([-0.5,0.5])/J
        M+=wi*J*np.outer(Nv,Nv); K+=wi*J*np.outer(dN,dN)
        Vg=float(np.dot(Nv,Vn)); Vi+=wi*J*Vg*np.outer(Nv,Nv)
    return K,M,Vi
def resolve(L,n):
    N=n+1; x=np.linspace(-L/2,L/2,N)
    K=np.zeros((N,N));M=np.zeros((N,N));V=np.zeros((N,N))
    for i in range(n):
        k,m,vi=mats1d(x[i],x[i+1])
        K[i:i+2,i:i+2]+=k; M[i:i+2,i:i+2]+=m; V[i:i+2,i:i+2]+=vi
    Kk=np.kron(K,M)+np.kron(M,K); Mm=np.kron(M,M)
    Vv=ALPHA*(np.kron(V,M)+np.kron(M,V))
    A=Kk+Vv
    livre=np.ones((N,N),bool); livre[0]=livre[-1]=False; livre[:,0]=livre[:,-1]=False
    f=np.where(livre.ravel())[0]
    Af=A[np.ix_(f,f)]; Mf=Mm[np.ix_(f,f)]
    Lc=np.linalg.cholesky(Mf)
    C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
    return np.linalg.eigvalsh(C)[0]
print("2D harmonico, caixa L=20, alpha=0.25.  Analitico ILIMITADO = 1.0")
print("   n     h        E1            erro      ordem")
ant=None
for n in [10,20,30,40,60,80]:
    e=resolve(20.0,n); err=abs(e-1.0); h=20.0/n
    o=(np.log(ant[1]/err)/np.log(ant[0]/h)) if ant else None
    print("  %3d  %.4f  %.10f  %.3e  %s"%(n,h,e,err,"%.2f"%o if o else "  -"))
    ant=(h,err)
print()
print("  artigo (n=40): 1.028585384356")
