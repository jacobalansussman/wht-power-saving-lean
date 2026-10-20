# Skeptic check (pure Python, exact rationals/Gaussian rationals, seconds, < 50 MB).
# (1) Build Strassen's algorithm for n x n matrices (n = 2, 4) as a straight-line ARITHMETIC circuit
#     (gates: +, -, *, scalar*), count its gates C(n). Fix the left factor X to constants, then linearise
#     every gate to (value at Y = 0, degree-1 part) as in the four-step claim, building the LINEAR circuit
#     gate by gate and counting its gates. Check the linear circuit computes Y -> X Y exactly and has
#     <= 3 C(n) gates. Then do the same with the right factor fixed, and compose to get
#     (H (x) H) vec(Y) on N = n^2 points with <= 6 C(n) gates.
# (2) A NON-bilinear circuit for XY (adds y11*y11*x11 - x11*y11*y11 to an output, and a product of two
#     Y-dependent gates that cancels): linearisation still yields the right linear map.
# (3) Bluestein: DFT_N = D * (cyclic convolution of length L >= 2N-1 with a fixed chirp) * D, where the
#     convolution circuit (schoolbook, an arithmetic circuit in a and b) is linearised with a fixed to the
#     chirp. Exact check over Q(zeta) is heavy, so numerically to 1e-12 for N = 2..9.
from fractions import Fraction as Fr
import cmath, math, random
random.seed(7)

class Circ:
    def __init__(s): s.g=[]          # gates: ('in',name) | ('const',c) | ('+',a,b) | ('-',a,b) | ('*',a,b) | ('s',c,a)
    def add(s,*t): s.g.append(t); return len(s.g)-1
    def cost(s): return sum(1 for t in s.g if t[0] in '+-*s')

def strassen(C,X,Y,n):
    # X, Y: n x n lists of gate ids; returns n x n gate ids for XY
    if n==1: return [[C.add('*',X[0][0],Y[0][0])]]
    h=n//2
    def blk(M,i,j): return [r[j*h:(j+1)*h] for r in M[i*h:(i+1)*h]]
    def ad(A,B,op): return [[C.add(op,A[i][j],B[i][j]) for j in range(h)] for i in range(h)]
    A11,A12,A21,A22=blk(X,0,0),blk(X,0,1),blk(X,1,0),blk(X,1,1)
    B11,B12,B21,B22=blk(Y,0,0),blk(Y,0,1),blk(Y,1,0),blk(Y,1,1)
    M1=strassen(C,ad(A11,A22,'+'),ad(B11,B22,'+'),h)
    M2=strassen(C,ad(A21,A22,'+'),B11,h)
    M3=strassen(C,A11,ad(B12,B22,'-'),h)
    M4=strassen(C,A22,ad(B21,B11,'-'),h)
    M5=strassen(C,ad(A11,A12,'+'),B22,h)
    M6=strassen(C,ad(A21,A11,'-'),ad(B11,B12,'+'),h)
    M7=strassen(C,ad(A12,A22,'-'),ad(B21,B22,'+'),h)
    C11=ad(ad(ad(M1,M4,'+'),M5,'-'),M7,'+'); C12=ad(M3,M5,'+'); C21=ad(M2,M4,'+'); C22=ad(ad(ad(M1,M2,'-'),M3,'+'),M6,'+')
    return [C11[i]+C12[i] for i in range(h)]+[C21[i]+C22[i] for i in range(h)]

def linearise(C,inputs_const,yvars):
    # Propagate (c0, lin) where lin is a dict var->coef; count linear gates as built:
    # +/-: 1 gate if both lin nonempty (else 0 or a free sign); s: 1 gate; *: u0*v1 + v0*u1 -> up to 3 gates.
    val={}; gates=0
    for i,t in enumerate(C.g):
        if t[0]=='in':
            nm=t[1]
            if nm in inputs_const: val[i]=(inputs_const[nm],{})
            else: val[i]=(0,{nm:1})
        elif t[0]=='const': val[i]=(t[1],{})
        elif t[0] in '+-':
            (a0,a1),(b0,b1)=val[t[1]],val[t[2]]; sg=1 if t[0]=='+' else -1
            d=dict(a1)
            for k,v in b1.items(): d[k]=d.get(k,0)+sg*v
            d={k:v for k,v in d.items() if v!=0}
            if a1 and b1: gates+=1
            elif b1 and sg==-1: gates+=1    # negation charged
            val[i]=(a0+sg*b0,d)
        elif t[0]=='s':
            c=t[1]; (a0,a1)=val[t[2]]
            if a1 and c not in (1,): gates+=1
            val[i]=(c*a0,{k:c*v for k,v in a1.items()})
        elif t[0]=='*':
            (a0,a1),(b0,b1)=val[t[1]],val[t[2]]
            d={}
            terms=0
            if b1 and a0!=0:
                terms+=1
                for k,v in b1.items(): d[k]=d.get(k,0)+a0*v
            if a1 and b0!=0:
                terms+=1
                for k,v in a1.items(): d[k]=d.get(k,0)+b0*v
            d={k:v for k,v in d.items() if v!=0}
            gates+= {0:0,1:1,2:3}[terms]
            val[i]=(a0*b0,d)
    return val,gates

def Hk(k):
    M=[[1]]
    for _ in range(k): M=[[M[i%len(M)][j%len(M)]*(1 if (i//len(M))*(j//len(M))==0 else -1) for j in range(2*len(M))] for i in range(2*len(M))]
    return M

out=[]
for k in (1,2):
    n=2**k; H=Hk(k)
    # left product: X fixed = H, Y variables
    C=Circ()
    X=[[C.add('in','x%d_%d'%(i,j)) for j in range(n)] for i in range(n)]
    Y=[[C.add('in','y%d_%d'%(i,j)) for j in range(n)] for i in range(n)]
    Z=strassen(C,X,Y,n); Cn=C.cost()
    val,gL=linearise(C,{'x%d_%d'%(i,j):Fr(H[i][j]) for i in range(n) for j in range(n)},None)
    okL=all(val[Z[i][j]][0]==0 and val[Z[i][j]][1]=={('y%d_%d'%(t,j)):H[i][t] for t in range(n) if H[i][t]!=0} for i in range(n) for j in range(n))
    # right product: Y fixed = H^T, X variables
    val2,gR=linearise(C,{'y%d_%d'%(i,j):Fr(H[j][i]) for i in range(n) for j in range(n)},None)
    okR=all(val2[Z[i][j]][0]==0 and val2[Z[i][j]][1]=={('x%d_%d'%(i,t)):H[j][t] for t in range(n) if H[j][t]!=0} for i in range(n) for j in range(n))
    # compose: (H Y) H^T = vec -> (H (x) H) vec Y ; check on random integer Y exactly
    Yv=[[random.randint(-5,5) for _ in range(n)] for _ in range(n)]
    HY=[[sum(val[Z[i][j]][1].get('y%d_%d'%(t,j),0)*Yv[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    R=[[sum(val2[Z[i][j]][1].get('x%d_%d'%(i,t),0)*HY[i][t] for t in range(n)) for j in range(n)] for i in range(n)]
    HH=[[H[a//n][b//n]*H[a%n][b%n] for b in range(n*n)] for a in range(n*n)]
    vecY=[Yv[i][j] for i in range(n) for j in range(n)]
    ref=[sum(HH[a][b]*vecY[b] for b in range(n*n)) for a in range(n*n)]
    okC=[R[i][j] for i in range(n) for j in range(n)]==ref
    out.append("Strassen n=%d: C(n)=%d gates; left-fixed linearised: correct=%s, %d gates; right-fixed: correct=%s, %d gates; total %d <= 6C(n)=%d: %s; (H(x)H)vecY exact: %s"
               %(n,Cn,okL,gL,okR,gR,gL+gR,6*Cn,gL+gR<=6*Cn,okC))

# (2) non-bilinear circuit: schoolbook XY for n=2 plus junk of degree 3 and a cancelling product of Y-gates
C=Circ(); n=2
X=[[C.add('in','x%d_%d'%(i,j)) for j in range(n)] for i in range(n)]
Y=[[C.add('in','y%d_%d'%(i,j)) for j in range(n)] for i in range(n)]
Z=[[None]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        acc=C.add('*',X[i][0],Y[0][j])
        for t in range(1,n): acc=C.add('+',acc,C.add('*',X[i][t],Y[t][j]))
        Z[i][j]=acc
one=C.add('const',Fr(1))
p=C.add('+',Y[0][0],one); q=C.add('-',Y[0][0],one)       # (y+1)(y-1) - y*y + 1 = 0
junk=C.add('+',C.add('-',C.add('*',p,q),C.add('*',Y[0][0],Y[0][0])),one)
junk2=C.add('-',C.add('*',C.add('*',Y[0][0],Y[0][0]),X[0][0]),C.add('*',X[0][0],C.add('*',Y[0][0],Y[0][0])))
Z[0][0]=C.add('+',Z[0][0],C.add('+',junk,junk2))
A=[[Fr(3),Fr(-1)],[Fr(2),Fr(5)]]
val,g=linearise(C,{'x%d_%d'%(i,j):A[i][j] for i in range(n) for j in range(n)},None)
ok=all(val[Z[i][j]][0]==0 and val[Z[i][j]][1]=={('y%d_%d'%(t,j)):A[i][t] for t in range(n) if A[i][t]!=0} for i in range(n) for j in range(n))
out.append("non-bilinear circuit (degree-3 junk + cancelling (y+1)(y-1)-y^2+1): C=%d gates, linearised correct=%s, %d gates <= 3C=%d: %s"%(C.cost(),ok,g,3*C.cost(),g<=3*C.cost()))

# (3) Bluestein through a linearised convolution circuit (numerical)
def conv_circuit(L):
    C=Circ(); a=[C.add('in','a%d'%i) for i in range(L)]; b=[C.add('in','b%d'%i) for i in range(L)]
    outg=[]
    for k in range(L):
        acc=None
        for i in range(L):
            t=C.add('*',a[i],b[(k-i)%L]); acc=t if acc is None else C.add('+',acc,t)
        outg.append(acc)
    return C,outg
for N in range(2,10):
    L=2*N-1
    w=lambda e: cmath.exp(-1j*math.pi*e/N)        # w^(e) with w = exp(-i pi/N), so w^(j^2) = exp(-i pi j^2/N)
    chirp=[0]*L
    for d in range(-(N-1),N): chirp[d%L]=w(-d*d)   # c_d = conj(w^{d^2})
    C,outg=conv_circuit(L)
    val,g=linearise(C,{'a%d'%i:chirp[i] for i in range(L)},None)
    x=[complex(random.gauss(0,1),random.gauss(0,1)) for _ in range(N)]
    bin_=[x[i]*w(i*i) if i<N else 0 for i in range(L)]
    conv=[sum(val[outg[k]][1].get('b%d'%i,0)*bin_[i] for i in range(L)) for k in range(L)]
    X=[w(k*k)*conv[k] for k in range(N)]
    ref=[sum(x[t]*cmath.exp(-2j*math.pi*t*k/N) for t in range(N)) for k in range(N)]
    err=max(abs(u-v) for u,v in zip(X,ref))
    out.append("Bluestein N=%d via linearised cyclic convolution of length %d: max err %.1e; linear gates %d <= 3*C_conv=%d: %s"%(N,L,err,g,3*C.cost(),g<=3*C.cost()))
print("\n".join(out))
