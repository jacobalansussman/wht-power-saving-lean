# Check 2 (exact, pure Python): the L0 log-exponent is the same for every nonmonomial 2x2 kernel.
# Words: (i) H = U * diag(1,-2) * P U P            (2 calls to the zeta/shear kernel U)
#        (ii) H D1 H D2 H = lower triangular with all three entries nonzero, D1 = diag(1,2), D2 = diag(1,-3),
#             hence U = monomial * (H D1 H D2 H) * monomial   (3 calls to H)
# and by the mixed-product rule (XY)^{(x)k} = X^{(x)k} Y^{(x)k}, the same words hold for k-th Kronecker powers,
# with monomial factors costing <= 2^k scalar gates each. Verified exactly for k = 1..4.
from fractions import Fraction as Fr
def mm(A,B): return [[sum(A[i][t]*B[t][j] for t in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def kron(A,B):
    return [[A[i//len(B)][j//len(B[0])]*B[i%len(B)][j%len(B[0])] for j in range(len(A[0])*len(B[0]))] for i in range(len(A)*len(B))]
def kp(A,k):
    M=[[Fr(1)]]
    for _ in range(k): M=kron(M,A)
    return M
def prod(*Ms):
    R=Ms[0]
    for M in Ms[1:]: R=mm(R,M)
    return R
F=lambda rows:[[Fr(x) for x in r] for r in rows]
H=F([[1,1],[1,-1]]); U=F([[1,0],[1,1]]); P=F([[0,1],[1,0]])
D=F([[1,0],[0,-2]]); D1=F([[1,0],[0,2]]); D2=F([[1,0],[0,-3]])
W=prod(H,D1,H,D2,H); print("H D1 H D2 H =", [[str(x) for x in r] for r in W])
# W = [[6,0],[-10,8]] = diag(6,8) * [[1,0],[-10/8,1]] ; and [[1,0],[c,1]] = diag(1,c) U diag(1,1/c)
a,c=W[0][0],W[1][0]; d=W[1][1]
c0=c/d
Dl=[[1/a,Fr(0)],[Fr(0),1/(d*c0)]]; Dr=[[Fr(1),Fr(0)],[Fr(0),c0]]
print("U == Dl W Dr :", prod(Dl,W,Dr)==U)
print("H == U D P U P :", prod(U,D,P,U,P)==H)
for k in range(1,5):
    ok1 = kp(H,k)==prod(kp(U,k),kp(D,k),kp(P,k),kp(U,k),kp(P,k))
    ok2 = kp(U,k)==prod(kp(Dl,k),kp(H,k),kp(D1,k),kp(H,k),kp(D2,k),kp(H,k),kp(Dr,k))
    print("k=%d  H^k via 2 calls to U^k: %s   U^k via 3 calls to H^k: %s" % (k,ok1,ok2))
print("=> L(H^k) <= 2 L(U^k) + 3*2^k and L(U^k) <= 3 L(H^k) + 4*2^k: WHT, Yates zeta transform and disjointness (= P U) share one log-exponent at L0.")
