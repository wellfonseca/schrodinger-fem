import numpy as np
ALPHA=0.25
def mats1d(a,b):
    """K, M (Gauss-2) e V_int, V_exa para o elemento [a,b]. V=x^2."""
    h=b-a; J=h/2
    g2=np.array([-1/np.sqrt(3),1/np.sqrt(3)]); w2=np.array([1.0,1.0])
    g3=np.array([-np.sqrt(3/5),0,np.sqrt(3/5)]); w3=np.array([5/9,8/9,5/9])
    K=np.zeros((2,2));M=np.zeros((2,2));Vi=np.zeros((2,2));Ve=np.zeros((2,2))
    Vn=np.array([a*a,b*b])
    for xi,wi in zip(g2,w2):
        Nv=np.array([(1-xi)/2,(1+xi)/2]); dN=np.array([-0.5,0.5])/J
        M+=wi*J*np.outer(Nv,Nv); K+=wi*J*np.outer(dN,dN)
        Vg=float(np.dot(Nv,Vn))
        Vi+=wi*J*Vg*np.outer(Nv,Nv)          # V interpolado (como Elmer)
    for xi,wi in zip(g3,w3):
        Nv=np.array([(1-xi)/2,(1+xi)/2]); xx=(a+b)/2+J*xi
        Ve+=wi*J*xx*xx*np.outer(Nv,Nv)       # V exato (Gauss-3)
    return K,M,Vi,Ve
def resolve(L,n,Vtipo):
    N=n+1; x=np.linspace(-L/2,L/2,N)
    K=np.zeros((N,N));M=np.zeros((N,N));V=np.zeros((N,N))
    for i in range(n):
        k,m,vi,ve=mats1d(x[i],x[i+1])
        K[i:i+2,i:i+2]+=k; M[i:i+2,i:i+2]+=m
        V[i:i+2,i:i+2]+= (vi if Vtipo=='int' else ve)
    K1,M1,V1=K,M,V
    I=np.eye(N)
    Kk=(np.kron(np.kron(K1,M1),M1)+np.kron(np.kron(M1,K1),M1)+np.kron(np.kron(M1,M1),K1))
    Mm=np.kron(np.kron(M1,M1),M1)
    Vv=ALPHA*(np.kron(np.kron(V1,M1),M1)+np.kron(np.kron(M1,V1),M1)+np.kron(np.kron(M1,M1),V1))
    A=Kk+Vv
    livre=np.ones((N,N,N),bool)
    livre[0]=livre[-1]=False; livre[:,0]=livre[:,-1]=False; livre[:,:,0]=livre[:,:,-1]=False
    f=np.where(livre.ravel())[0]
    Af=A[np.ix_(f,f)]; Mf=Mm[np.ix_(f,f)]
    Lc=np.linalg.cholesky(Mf)
    C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
    return np.linalg.eigvalsh(C)[0]
print("3D harmonico, caixa L=12, alpha=0.25.  Analitico ILIMITADO = 1.5")
print("Se a sequencia converge para 1.5 -> o erro e' DISCRETIZACAO, nao truncamento")
print()
for tipo,rot in [('int','V interpolado (como Elmer)'),('exa','V exato (Gauss-3)')]:
    print("  --- %s ---" % rot)
    print("     n     h        E1            erro      ordem")
    ant=None
    for n in [6,10,14,18,20]:
        e=resolve(12.0,n,tipo); err=abs(e-1.5)/1.5; h=12.0/n
        o=(np.log(ant[1]/err)/np.log(ant[0]/h)) if ant else None
        print("    %3d  %.4f  %.10f  %.3e  %s"%(n,h,e,err,"%.2f"%o if o else "  -"))
        ant=(h,err)
    print()
