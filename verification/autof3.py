import numpy as np
g3=np.array([-0.7745966692414834,0.0,0.7745966692414834]); w3=np.array([5.0/9,8.0/9,5.0/9])
def mede(L,n):
    N=n+1; x=np.linspace(-L/2,L/2,N); h=L/n
    K=np.zeros((N,N));M=np.zeros((N,N))
    for i in range(n):
        K[i:i+2,i:i+2]+=np.array([[1,-1],[-1,1]])/h
        M[i:i+2,i:i+2]+=h/6*np.array([[2,1],[1,2]])
    f=np.arange(1,N-1); Af=K[np.ix_(f,f)]; Mf=M[np.ix_(f,f)]
    Lc=np.linalg.cholesky(Mf); C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
    val,y=np.linalg.eigh(C); u=np.linalg.solve(Lc.T,y[:,0]); E=val[0]
    u=u/np.sqrt(u@(Mf@u))
    psi=np.zeros(N); psi[f]=u
    psia=lambda xx: np.sqrt(2.0/L)*np.cos(np.pi*xx/L)
    dpsia=lambda xx:-np.sqrt(2.0/L)*(np.pi/L)*np.sin(np.pi*xx/L)
    # SIGN ALIGNMENT: psi(0) must be positive
    ic=len(x)//2
    if psi[ic] < 0: psi=-psi
    L2=0.0; H1s=0.0
    for i in range(n):
        a,b=x[i],x[i+1]; J=(b-a)/2
        for xi,wi in zip(g3,w3):
            xx=(a+b)/2+J*xi; Nv=np.array([(1-xi)/2,(1+xi)/2])
            uh=float(Nv@np.array([psi[i],psi[i+1]])); duh=(psi[i+1]-psi[i])/h
            L2+=wi*J*(uh-psia(xx))**2; H1s+=wi*J*(duh-dpsia(xx))**2
    return E, np.sqrt(L2), np.sqrt(L2+H1s)
print("  1D infinite well L=40 -- EIGENFUNCTION error vs exact function (sign aligned)")
print("     n     h        ||e||_L2     order    ||e||_H1     order")
ant=None
for n in [25,50,100,200,400]:
    E,L2,H1=mede(40.0,n); h=40.0/n
    if ant:
        print("  %4d  %.4f  %.4e   %.2f     %.4e   %.2f" % (n,h,L2,
            np.log(ant[1]/L2)/np.log(ant[0]/h), H1, np.log(ant[2]/H1)/np.log(ant[0]/h)))
    else:
        print("  %4d  %.4f  %.4e     -       %.4e     -" % (n,h,L2,H1))
    ant=(h,L2,H1)
