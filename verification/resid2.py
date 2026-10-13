import numpy as np
g2=np.array([-1/np.sqrt(3),1/np.sqrt(3)]); w2=np.array([1.0,1.0])
def m1d(a,b):
    """K, M e V1 = int Ni x^2 Nj  (potencial x^2, sem alpha)."""
    h=b-a; J=h/2; K=np.zeros((2,2));M=np.zeros((2,2));V=np.zeros((2,2))
    Vn=np.array([a*a,b*b])
    for xi,wi in zip(g2,w2):
        Nv=np.array([(1-xi)/2,(1+xi)/2]); dN=np.array([-0.5,0.5])/J
        M+=wi*J*np.outer(Nv,Nv); K+=wi*J*np.outer(dN,dN)
        V+=wi*J*float(np.dot(Nv,Vn))*np.outer(Nv,Nv)
    return K,M,V
def blk(ms,N):
    A=np.zeros((N,N))
    for i,m in enumerate(ms): A[i:i+2,i:i+2]+=m
    return A
def caso(L,n,dim,alpha):
    N=n+1; x=np.linspace(-L/2,L/2,N)
    Ks=[];Ms=[];Vs=[]
    for i in range(n):
        k,m,v=m1d(x[i],x[i+1]); Ks.append(k);Ms.append(m);Vs.append(v)
    K1,M1,V1=blk(Ks,N),blk(Ms,N),blk(Vs,N)
    if dim==1:
        return K1+alpha*V1, M1, np.ones(N,bool)
    if dim==2:
        K=np.kron(K1,M1)+np.kron(M1,K1); M=np.kron(M1,M1)
        V=alpha*(np.kron(V1,M1)+np.kron(M1,V1))
        livre=np.ones((N,N),bool); livre[0]=livre[-1]=False
        livre[:,0]=livre[:,-1]=False
        return K+V, M, livre
    K=(np.kron(np.kron(K1,M1),M1)+np.kron(np.kron(M1,K1),M1)+np.kron(np.kron(M1,M1),K1))
    M=np.kron(np.kron(M1,M1),M1)
    V=alpha*(np.kron(np.kron(V1,M1),M1)+np.kron(np.kron(M1,V1),M1)+np.kron(np.kron(M1,M1),V1))
    livre=np.ones((N,N,N),bool); livre[0]=livre[-1]=False
    livre[:,0]=livre[:,-1]=False; livre[:,:,0]=livre[:,:,-1]=False
    return K+V, M, livre
for rot,(L,n,dim,al,el) in [
   ("1D poco infinito",(40.0,250,1,0.0,6.168583925347e-03)),
   ("2D harmonico",    (20.0,40,2,0.25,1.028585384356e+00)),
   ("3D harmonico",    (12.0,14,3,0.25,1.625507320678e+00))]:
    A,M,livre=caso(L,n,dim,al)
    livre[0]=livre[-1]=False
    if dim>=2: livre[:,0]=livre[:,-1]=False
    if dim==3: livre[:,:,0]=livre[:,:,-1]=False
    f=np.where(livre.ravel())[0]
    print("  %-18s A=%s  f=%d" % (rot, A.shape, len(f)), flush=True)
    Af=A[np.ix_(f,f)]; Mf=M[np.ix_(f,f)]
    Lc=np.linalg.cholesky(Mf); C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
    val,y=np.linalg.eigh(C); u=np.linalg.solve(Lc.T,y[:,0]); E=val[0]
    u=u/np.sqrt(u@(Mf@u))
    r=Af@u - E*(Mf@u)
    print("      E=%.12f  |E-Elmer|=%.1e  res.norm=%.2e  |orto-1|=%.1e"
          % (E, abs(E-el), np.linalg.norm(r)/np.linalg.norm(Af@u), abs(u@(Mf@u)-1)), flush=True)
