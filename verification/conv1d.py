import numpy as np
def mats1d(a,b):
    h=b-a; g=np.array([-np.sqrt(3/5),0,np.sqrt(3/5)]); w=np.array([5/9,8/9,5/9])
    K=np.zeros((2,2));M=np.zeros((2,2));V=np.zeros((2,2))
    for xi,wi in zip(g,w):
        N=np.array([(1-xi)/2,(1+xi)/2]); x=(a+b)/2+h/2*xi; J=h/2
        M+=wi*J*np.outer(N,N); V+=wi*J*x*x*np.outer(N,N)
        dN=np.array([-0.5,0.5])/J; K+=wi*J*np.outer(dN,dN)
    return K,M,V
def E1(L,n,alpha):
    N=n+1; x=np.linspace(-L/2,L/2,N)
    K=np.zeros((N,N));M=np.zeros((N,N));V=np.zeros((N,N))
    for i in range(n):
        k,m,v=mats1d(x[i],x[i+1])
        K[i:i+2,i:i+2]+=k; M[i:i+2,i:i+2]+=m; V[i:i+2,i:i+2]+=alpha*v
    A=K+V
    f=np.arange(1,N-1)
    Af=A[np.ix_(f,f)]; Mf=M[np.ix_(f,f)]
    Lc=np.linalg.cholesky(Mf)
    C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
    return np.linalg.eigvalsh(C)[0]
print("1D harmonic:  V = x^2  (alpha=1),  box L=20  -> analytic E0 = 1.0")
print("   n     h        E1            error       order")
ant=None
for n in [25,50,100,200,400,800]:
    e=E1(20.0,n,1.0); err=e-1.0; h=20.0/n
    o=(np.log(ant[1]/err)/np.log(ant[0]/h)) if ant else None
    print("  %4d  %.4f  %.10f  %.3e  %s"%(n,h,e,err,"%.2f"%o if o else "  -"))
    ant=(h,err)
print()
print("  paper (Elmer, 400 elem): 1.0005728920")
print("  mine (400 elem)         : %.10f" % E1(20.0,400,1.0))
