"""Skeptic toy: restriction of scalars turns a unit with coefficients in L = Q(t), t^3 = 2 (non-abelian
Galois closure) into a unit with RATIONAL mixing, same blocks per live slot.  Engine-style semantics:
Data = roles -> (F_2)^f -> C; a fixed letter M (2x2 over Q(i)) on role l acts as M^{(x)f}; localM g mixes
roles pointwise.  Exact arithmetic over Q(i) for the restricted word; floats for the original via embeddings."""
from fractions import Fraction as Fr
import itertools, random, cmath
random.seed(7)
class Gi:                      # Q(i), exact
    __slots__=('a','b')
    def __init__(s,a,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o): o=o if isinstance(o,Gi) else Gi(o); return Gi(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __sub__(s,o): o=o if isinstance(o,Gi) else Gi(o); return Gi(s.a-o.a,s.b-o.b)
    def __mul__(s,o): o=o if isinstance(o,Gi) else Gi(o); return Gi(s.a*o.a-s.b*o.b,s.a*o.b+s.b*o.a)
    __rmul__=__mul__
    def __eq__(s,o): o=o if isinstance(o,Gi) else Gi(o); return s.a==o.a and s.b==o.b
    def c(s): return complex(float(s.a),float(s.b))
I=Gi(0,1); h=Fr(1,2)
ca=Gi(h,h); cb=Gi(h,-h)
K=[[ca,cb],[cb,ca]]                 # upstream littleC
P=[[Gi(1),Gi(0)],[Gi(0),I]]         # phaseMat z=1 c=0: diag(1,i)
# L = Q(t), t^3=2 ; element = (a0,a1,a2)
def Lmul(x,y):
    c=[Fr(0)]*5
    for i in range(3):
        for j in range(3): c[i+j]+=x[i]*y[j]
    return (c[0]+2*c[3], c[1]+2*c[4], c[2])
def Linv(x):                        # solve x*y=1 by the regular matrix
    R=[[Lmul(x,e)[r] for e in ((1,0,0),(0,1,0),(0,0,1))] for r in range(3)]
    # Gaussian elimination over Q
    A=[R[r][:]+[Fr(1 if r==0 else 0)] for r in range(3)]
    for c in range(3):
        p=next(r for r in range(c,3) if A[r][c]!=0); A[c],A[p]=A[p],A[c]
        A[c]=[v/A[c][c] for v in A[c]]
        for r in range(3):
            if r!=c and A[r][c]!=0: A[r]=[a-A[r][c]*b for a,b in zip(A[r],A[c])]
    return tuple(A[r][3] for r in range(3))
def Ladd(x,y): return tuple(a+b for a,b in zip(x,y))
def Lsub(x,y): return tuple(a-b for a,b in zip(x,y))
def rnd(): return tuple(Fr(random.randint(-5,5),random.randint(1,4)) for _ in range(3))
# ---- an L-unit: roles {0,1,2}, live {0}; word = localM g1 ; K@0 ; P@1 ; K@2 ; localM g2
# output role 0 = sum_j g2[0][j] M_j g1[j][0] x ; need p0+p2 = 1 (K class), p1 = 0 (P class), p_j = g2[0][j]*g1[j][0]
g1=[[rnd() for _ in range(3)] for _ in range(3)]
g1[1][0]=(Fr(0),)*3
g2=[[rnd() for _ in range(3)] for _ in range(3)]
g2[0][2]=Lmul(Lsub((Fr(1),Fr(0),Fr(0)),Lmul(g2[0][0],g1[0][0])),Linv(g1[2][0]))
print('irrational entries (t-coefficients nonzero):', sum(1 for g in (g1,g2) for row in g for x in row if x[1]!=0 or x[2]!=0), 'of 18')
mats={0:K,1:P,2:K}
def kron_pow(M,f):
    out=[[Gi(1)]]
    for _ in range(f):
        out=[[out[i][j]*M[a][b] for j in range(len(out)) for b in range(2)] for i in range(len(out)) for a in range(2)]
    return out
def mv(M,v): return [sum((M[i][j]*v[j] for j in range(len(v))),Gi(0)) for i in range(len(M))]
# ---- restricted (Weil) word: roles (r,b), rational mixing = regular representation
basis=((Fr(1),Fr(0),Fr(0)),(Fr(0),Fr(1),Fr(0)),(Fr(0),Fr(0),Fr(1)))
def reg(g):   # g: 3x3 over L -> 9x9 over Q
    G=[[Fr(0)]*9 for _ in range(9)]
    for r in range(3):
        for j in range(3):
            for bp in range(3):
                prod=Lmul(g[r][j],basis[bp])
                for b in range(3): G[3*r+b][3*j+bp]=prod[b]
    return G
G1,G2=reg(g1),reg(g2)
def run_restricted(f,v):            # v: dict role(0..8) -> list of 2^f Gi
    n=2**f
    def mix(G,x): return {R:[sum((Gi(G[R][S])*x[S][i] for S in range(9) if G[R][S]!=0),Gi(0)) for i in range(n)] for R in range(9)}
    x=mix(G1,v)
    for R in range(9): x[R]=mv(kron_pow(mats[R//3],f),x[R])
    return mix(G2,x)
ok=True
for f in (1,2,3):
    n=2**f
    for trial in range(3):
        v={R:([Gi(Fr(random.randint(-9,9),random.randint(1,5)),Fr(random.randint(-9,9),random.randint(1,5))) for _ in range(n)] if R//3==0 else [Gi(0)]*n) for R in range(9)}
        out=run_restricted(f,v)
        Kf=kron_pow(K,f)
        for b in range(3):
            if not all(a==c for a,c in zip(out[b],mv(Kf,v[b]))): ok=False
print('restricted RATIONAL-mixing word: live outputs (roles (0,b), b=0,1,2) equal K^(x)f input, exactly, f=1,2,3:',ok)
print('mixing entries all rational:', all(isinstance(x,Fr) for G in (G1,G2) for row in G for x in row))
print('paid blocks / live slots: original 2/1 ; restricted', 2*3, '/', 3)
# ---- original word under the three embeddings of L (t -> 2^(1/3) w^k): each must satisfy the system
for k in range(3):
    tt=2**(1/3)*cmath.exp(2j*cmath.pi*k/3)
    ev=lambda x: complex(x[0])+complex(x[1])*tt+complex(x[2])*tt*tt
    pK=ev(Lmul(g2[0][0],g1[0][0]))+ev(Lmul(g2[0][2],g1[2][0])); pP=ev(Lmul(g2[0][1],g1[1][0]))
    print(f'embedding {k}: K-coefficient {pK:.3g} (target 1), P-coefficient {abs(pP):.1e} (target 0)')
