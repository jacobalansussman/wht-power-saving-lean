# Check 1 (pure Python, no numpy): the MM -> transform reduction (four-step), on small cases.
# (a) WHT on 4^k entries = two products with the fixed 2^k x 2^k matrix H_k:
#     (H_k (x) H_k) vec(X) = vec(H_k X H_k^T)   (row-major vec), exact integer arithmetic.
# (b) DFT on n^2 entries = two products with F_n plus n^2 twiddles (Gentleman-Sande / Bailey four-step).
# (c) Gate count of Strassen's algorithm with the left factor fixed, linearised (scalar circuit for Y -> A Y):
#     shows the MM route beats N log N only if MM costs n^2 (log n)^c with c < 1.
import cmath, math, random
random.seed(1)
def matmul(A,B):
    n,m,p = len(A),len(B),len(B[0])
    return [[sum(A[i][t]*B[t][j] for t in range(m)) for j in range(p)] for i in range(n)]
def kron(A,B):
    return [[A[i//len(B)][j//len(B[0])]*B[i%len(B)][j%len(B[0])] for j in range(len(A[0])*len(B[0]))] for i in range(len(A)*len(B))]
def T(A): return [list(r) for r in zip(*A)]
def Hk(k):
    M=[[1]]
    for _ in range(k): M = kron(M,[[1,1],[1,-1]])
    return M
for k in range(1,4):
    n=2**k; X=[[random.randint(-9,9) for _ in range(n)] for _ in range(n)]
    vec=[X[i][j] for i in range(n) for j in range(n)]
    lhs=[sum(r[t]*vec[t] for t in range(n*n)) for r in kron(Hk(k),Hk(k))]
    R=matmul(matmul(Hk(k),X),T(Hk(k))); rhs=[R[i][j] for i in range(n) for j in range(n)]
    print("WHT four-step k=%d N=%d exact equal: %s" % (k,n*n,lhs==rhs))
def F(n): return [[cmath.exp(-2j*math.pi*a*b/n) for b in range(n)] for a in range(n)]
for n in [2,3,4,5,6]:
    N=n*n; x=[complex(random.gauss(0,1),random.gauss(0,1)) for _ in range(N)]
    X=[[x[n*a+b] for b in range(n)] for a in range(n)]        # input t = n a + b
    Y=matmul(F(n),X)                                          # Y[c][b] = sum_a w_n^{ac} X[a][b]
    Yt=[[Y[c][b]*cmath.exp(-2j*math.pi*c*b/N) for b in range(n)] for c in range(n)]
    Z=matmul(Yt,F(n))                                         # Z[c][d]
    out=[Z[c][d] for d in range(n) for c in range(n)]         # output f = c + n d
    ref=[sum(F(N)[f][t]*x[t] for t in range(N)) for f in range(N)]
    print("DFT four-step n=%d N=%d max err %.1e" % (n,N,max(abs(u-v) for u,v in zip(out,ref))))
# (c) Strassen with A fixed: of its 18 block additions, 5 are A-only (free, precomputed), 13 touch Y or outputs.
def lin(k): return 1 if k==0 else 7*lin(k-1)+13*4**(k-1)
for k in [2,4,6,8,10,12,16,20]:
    N=4**k; S=2*lin(k)
    print("k=%2d N=4^%-2d Strassen-route linear gates %.3e   N log2 N %.3e   ratio %.2f" % (k,k,S,N*math.log2(N),S/(N*math.log2(N))))
# (d) the condition on MM: if MM(n) <= C n^2 (log2 n)^c then WHT_N <= 6 C N (log2 N / 2)^c (linearisation factor 3, two products)
for c in [1.0,0.9,0.5]:
    print("MM in C n^2 (log n)^%.1f  ==>  WHT/DFT in O(N (log N)^%.1f), saving s >= %.1f" % (c,c,1-c))
