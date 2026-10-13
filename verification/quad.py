import numpy as np
def E1(L,n,alpha,ngauss):
    N=n+1; x=np.linspace(-L/2,L/2,N)
    if ngauss==1: g,w=np.array([0.0]),np.array([2.0])
    if ngauss==2: g,w=np.array([-1/np.sqrt(3),1/np.sqrt(3)]),np.array([1.0,1.0])
    if ngauss==3: g,w=np.array([-np.sqrt(3/5),0,np.sqrt(3/5)]),np.array([5/9,8/9,5/9])
    K=np.zeros((N,N));M=np.zeros((N,N));V=np.zeros((N,N))
    for i in range(n):
        a=x[i]; b=x[i+1]; h=b-a; J=h/2
        Kl=np.zeros((2,2));Ml=np.zeros((2,2));Vl=np.zeros((2,2))
        for xi,wi in zip(g,w):
            Nv=np.array([(1-xi)/2,(1+xi)/2]); xx=(a+b)/2+J*xi
            Ml+=wi*J*np.outer(Nv,Nv); Vl+=wi*J*xx*xx*np.outer(Nv,Nv)
            dN=np.array([-0.5,0.5])/J; Kl+=wi*J*np.outer(dN,dN)
        K[i:i+2,i:i+2]+=Kl; M[i:i+2,i:i+2]+=Ml; V[i:i+2,i:i+2]+=alpha*Vl
    A=K+V; f=np.arange(1,N-1)
    Af=A[np.ix_(f,f)]; Mf=M[np.ix_(f,f)]
    Lc=np.linalg.cholesky(Mf)
    C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
    return np.linalg.eigvalsh(C)[0]
print("1D harmonico, L=20, n=400, alpha=1  (analitico = 1.0)")
print("  artigo/Elmer : 1.000572892257")
print()
for ng in [1,2,3,5]:
    e=E1(20.0,400,1.0,ng)
    print("  Gauss %d ponto(s): %.12f   erro %.3e" % (ng,e,e-1.0))
