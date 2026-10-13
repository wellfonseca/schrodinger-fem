import numpy as np
g5=np.array([-0.9061798459386640,-0.5384693101056831,0.0,0.5384693101056831,0.9061798459386640])
w5=np.array([0.2369268850561891,0.4786286704993665,0.5688888888888889,0.4786286704993665,0.2369268850561891])
def calc(n):
    L=20.0; N=n+1; x=np.linspace(-L/2,L/2,N)
    K=np.zeros((N,N));M=np.zeros((N,N));V=np.zeros((N,N))
    for i in range(n):
        a,b=x[i],x[i+1]; h=b-a; J=h/2
        for xi,wi in zip(g5,w5):
            Nv=np.array([(1-xi)/2,(1+xi)/2]); dN=np.array([-0.5,0.5])/J
            K[i:i+2,i:i+2]+=wi*J*np.outer(dN,dN); M[i:i+2,i:i+2]+=wi*J*np.outer(Nv,Nv)
            Vg=float(np.dot(Nv,[a*a,b*b]))
            V[i:i+2,i:i+2]+=wi*J*Vg*np.outer(Nv,Nv)
    f=np.arange(1,N-1)
    Af=(K+V)[np.ix_(f,f)]; Mf=M[np.ix_(f,f)]; Kf=K[np.ix_(f,f)]; Vf=V[np.ix_(f,f)]
    Lc=np.linalg.cholesky(Mf)
    C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
    val,vec=np.linalg.eigh(C)
    y=vec[:,0]
    u=np.linalg.solve(Lc.T,y)          # <<< CORRECT BACK-TRANSFORMATION
    nrm=u@(Mf@u)
    E_rq=(u@(Af@u))/nrm
    T=(u@(Kf@u))/nrm; Vv=(u@(Vf@u))/nrm
    return val[0],E_rq,T,Vv
print("1D harmonic - virial with the FEM matrices (correct back-transformation)")
print("   n     E (eigh)       E_Rayleigh     <T>            <V>            <T>+<V>       ratio")
for n in [100,200,400,800,1600]:
    E,Erq,T,Vv=calc(n)
    print("  %4d  %.10f  %.10f  %.10f  %.10f  %.10f  %.8f"%(n,E,Erq,T,Vv,T+Vv,T/Vv))
print()
print("  exact: <T> = <V> = 0.5, ratio = 1")
