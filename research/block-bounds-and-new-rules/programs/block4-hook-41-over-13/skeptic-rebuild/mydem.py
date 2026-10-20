# Independent skeptic recomputation of dem = f(K) for the hook certificate N = span{e_x (x) e_T : x in T}
# in V (x) Lambda^k V, K = H^{(x)m} (Sylvester, +-1; equals sqrt(X)^{(x)m} up to diagonal phases, which keep N).
# dem = rank of A = P_{N^c} K |_N (integer matrix).  Decomposition (derived by me):
#  k even: dem = n * rank M_alpha (any alpha; all equal), M_alpha[S,T] = det H[S,T], rows S with alpha+sum(S) not in S, cols T with 0 in T.
#  k odd : dem = sum_{alpha,beta} rank M_{alpha,beta}, M_{ab}[S,T] = det H[S,T], rows 0 notin S, sum S = alpha; cols 0 in T, sum T = beta.
# Ranks mod a prime p (lower bound on rank over Q), optional random compression (rank(RA) <= rank A).
import sys  # needs numpy
import itertools, time, numpy as np
from math import comb
from functools import reduce
m=int(sys.argv[1]); k=int(sys.argv[2]); mode=sys.argv[3] if len(sys.argv)>3 else 'all'
n=2**m
H=np.array([[(-1)**bin(x&y).count('1') for y in range(n)] for x in range(n)],dtype=np.float64)
P=16777213  # prime < 2^24
def isprime(q):
    d=2
    while d*d<=q:
        if q%d==0: return False
        d+=1
    return True
assert isprime(P)
def xs(S): return reduce(lambda a,b:a^b,S,0)
def minors(Slist,Tlist):
    Sa=np.array(Slist,dtype=np.int64); out=np.zeros((len(Slist),len(Tlist)),dtype=np.int16)
    for j,T in enumerate(Tlist):
        sub=H[Sa][:,:,list(T)]  # (rows,k,k)
        d=np.linalg.det(sub); r=np.rint(d)
        assert np.max(np.abs(d-r))<1e-6
        out[:,j]=r.astype(np.int16)
    return out
def rank_modp(M,p=P):
    M=np.array(M,dtype=np.int64)%p
    r,c=M.shape; rk=0
    for col in range(c):
        if rk==r: break
        colv=M[rk:,col]; nz=np.nonzero(colv)[0]
        if len(nz)==0: continue
        piv=rk+nz[0]
        if piv!=rk: M[[rk,piv]]=M[[piv,rk]]
        inv=pow(int(M[rk,col]),p-2,p)
        M[rk,col:]=(M[rk,col:]*inv)%p
        f=M[rk+1:,col].copy(); nzf=np.nonzero(f)[0]
        if len(nzf):
            rows=rk+1+nzf
            M[rows,col:]=(M[rows,col:]-np.outer(f[nzf],M[rk,col:]))%p
        rk+=1
    return rk
def compress(A,p=P,seed=1,extra=8):
    r,c=A.shape
    if r<=c+extra: return A.astype(np.int64)
    rng=np.random.default_rng(seed); t=c+extra
    B=np.zeros((t,c),dtype=np.float64); ch=4096
    Bx=np.zeros((4,c),dtype=np.int64)   # exact int64 spot check of 4 compressed rows
    with np.errstate(all='ignore'):
        for s in range(0,r,ch):
            R=rng.integers(0,p,size=(t,min(ch,r-s)))
            Ac=A[s:s+ch]
            B=np.fmod(B+R.astype(np.float64)@Ac.astype(np.float64),p)
            Bx=(Bx+R[:4].astype(np.int64)@Ac.astype(np.int64))%p
    assert np.isfinite(B).all()
    B=np.mod(np.rint(B),p).astype(np.int64)
    assert (B[:4]==Bx).all(), "float compression not exact"
    return B
t0=time.time()
allS=list(itertools.combinations(range(n),k))
if k%2==0:
    alphas=[0] if mode=='all' else [int(a) for a in mode.split(',')]
    Tcols=[T for T in allS if 0 in T]
    for al in alphas:
        Srows=[S for S in allS if (al^xs(S)) not in S]
        A=minors(Srows,Tcols)
        rk=rank_modp(compress(A))
        print(f"m={m} k={k} alpha={al}: rows {len(Srows)} cols {len(Tcols)} rank_p {rk} -> n*rank = {n*rk}  ({time.time()-t0:.1f}s)",flush=True)
else:
    def blk(al,be):
        Srows=[S for S in allS if 0 not in S and xs(S)==al]
        Tcols=[T for T in allS if 0 in T and xs(T)==be]
        if not Srows or not Tcols: return 0,len(Srows),len(Tcols)
        A=minors(Srows,Tcols); return rank_modp(compress(A)),len(Srows),len(Tcols)
    if mode=='all':
        tot=0
        for al in range(n):
            for be in range(n):
                r,a,b=blk(al,be); tot+=r
        print(f"m={m} k={k}: ALL blocks, dem >= {tot}  ({time.time()-t0:.1f}s)")
    else:
        # 5 GL-orbit types: (0,0),(a,0),(0,b),(a,b: a.b=0),(a,b: a.b=1); plus random extra blocks to test the symmetry
        pc=lambda x: bin(x).count('1')
        reps={'00':(0,0),'a0':(1,0),'0b':(0,1),'ab0':(1,2),'ab1':(1,1)}
        size={'00':1,'a0':n-1,'0b':n-1,'ab0':(n-1)*(n//2-1),'ab1':(n-1)*(n//2)}
        tot=0
        for key,(al,be) in reps.items():
            r,a,b=blk(al,be); tot+=size[key]*r
            print(f"  type {key} (alpha={al},beta={be}) rows {a} cols {b} rank_p {r} x{size[key]}  ({time.time()-t0:.1f}s)",flush=True)
        if mode=='types+check':
            rng=np.random.default_rng(7)
            for _ in range(4):
                al,be=int(rng.integers(1,n)),int(rng.integers(1,n))
                key='ab1' if pc(al&be)%2 else 'ab0'
                r,a,b=blk(al,be); print(f"  check (alpha={al},beta={be}) type {key}: rank_p {r}",flush=True)
        print(f"m={m} k={k}: dem >= {tot} (by GL types)  ({time.time()-t0:.1f}s)")
print("delta upper bound", k*comb(n,k)-comb(n,k-1))
