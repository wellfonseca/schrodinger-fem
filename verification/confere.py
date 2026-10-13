import numpy as np
def mats1d(a,b,alpha=0.0):
    h=b-a; g=np.array([-np.sqrt(3/5),0,np.sqrt(3/5)]); w=np.array([5/9,8/9,5/9])
    K=np.zeros((2,2));M=np.zeros((2,2));V=np.zeros((2,2))
    for xi,wi in zip(g,w):
        N=np.array([(1-xi)/2,(1+xi)/2]); x=(a+b)/2+h/2*xi; J=h/2
        M+=wi*J*np.outer(N,N); V+=wi*J*x*x*np.outer(N,N)
        dN=np.array([-0.5,0.5])/J; K+=wi*J*np.outer(dN,dN)
    return K,M,V
def blk(mats,N):
    A=np.zeros((N,N))
    for i,m in enumerate(mats): A[i:i+2,i:i+2]+=m
    return A
def resolve(L,n,dim,alpha=0.25):
    N=n+1; x=np.linspace(-L/2,L/2,N)
    Ks=[];Ms=[];Vs=[]
    for i in range(n):
        k,m,v=mats1d(x[i],x[i+1]); Ks.append(k);Ms.append(m);Vs.append(v)
    K1,M1,V1=blk(Ks,N),blk(Ms,N),blk(Vs,N)
    if dim==1: K,M,V=K1,M1,V1
    if dim==2:
        K=np.kron(K1,M1)+np.kron(M1,K1); M=np.kron(M1,M1)
        V=alpha*(np.kron(V1,M1)+np.kron(M1,V1))
    if dim==3:
        K=(np.kron(np.kron(K1,M1),M1)+np.kron(np.kron(M1,K1),M1)+np.kron(np.kron(M1,M1),K1))
        M=np.kron(np.kron(M1,M1),M1)
        V=alpha*(np.kron(np.kron(V1,M1),M1)+np.kron(np.kron(M1,V1),M1)+np.kron(np.kron(M1,M1),V1))
    A=K+V
    if dim==1: livre=np.ones(N,bool); livre[0]=livre[-1]=False
    if dim==2:
        livre=np.ones((N,N),bool); livre[0]=livre[-1]=False; livre[:,0]=livre[:,-1]=False
    if dim==3:
        livre=np.ones((N,N,N),bool); livre[0]=livre[-1]=False
        livre[:,0]=livre[:,-1]=False; livre[:,:,0]=livre[:,:,-1]=False
    f=np.where(livre.ravel())[0]
    Af=A[np.ix_(f,f)]; Mf=M[np.ix_(f,f)]
    Lc=np.linalg.cholesky(Mf)
    C=np.linalg.solve(Lc,Af); C=np.linalg.solve(Lc,C.T).T; C=(C+C.T)/2
    val=np.linalg.eigvalsh(C)[:4]
    return val,len(f)
print("Case                  dim   L     n    E1(mine)       E1(paper)      dif")
casos=[("1D harmonic",1,20.0,400,1.000572892),
       ("2D harmonic",2,20.0,40, 1.028585384),
       ("3D harmonic",3,12.0,14, 1.625507321)]
for nome,d,L,n,art in casos:
    v,nd=resolve(L,n,d)
    print("  %-18s %d  %5.1f %4d  %.10f  %.10f  %+.2e"%(nome,d,L,n,v[0],art,v[0]-art))
print()
print("  analytic (infinite domain): 1D=1.0  2D=1.0  3D=1.5")
