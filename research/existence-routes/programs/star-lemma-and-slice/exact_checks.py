# Exact checks for the local hyperplane-rank problem (e3-slice6): ranks over Q with sympy.
from fractions import Fraction
def qrank(M):
    A=[[Fraction(v) for v in row] for row in M]; r=0
    rows=len(A); cols=len(A[0]) if rows else 0
    for c in range(cols):
        p=next((i for i in range(r,rows) if A[i][c]!=0),None)
        if p is None: continue
        A[r],A[p]=A[p],A[r]
        for i in range(rows):
            if i!=r and A[i][c]!=0:
                f=A[i][c]/A[r][c]; A[i]=[a-f*b for a,b in zip(A[i],A[r])]
        r+=1
    return r
def hyperplanes(n):
    N=1<<n; out=[]
    for h in range(1,N):
        for b in (0,1): out.append([x for x in range(N) if bin(h&x).count('1')%2==b])
    return out
def ranks(E, n):
    return [qrank([[E[i][j] for j in H] for i in H]) for H in hyperplanes(n)]
# (1) n=3, s=2: Singer permutation x -> alpha x in F_8 (perm_check.py's best: 0,2,4,6,3,1,7,5)
sig = (0, 2, 4, 6, 3, 1, 7, 5)
P = [[1 if sig[x] == y else 0 for y in range(8)] for x in range(8)]
r = ranks(P, 3); print("n=3 Singer permutation: det", "rank %d" % qrank(P), "max hyperplane rank", max(r), "ranks", r)
# (2) n=4, s=4: circulant E_xy = (-1)^phi(x+y), phi = x1x2 + x3x4 (bent)
phi = lambda x: ((x & 1) * (x >> 1 & 1) + (x >> 2 & 1) * (x >> 3 & 1)) % 2
C = [[(-1) ** phi(x ^ y) for y in range(16)] for x in range(16)]
r = ranks(C, 4); print("n=4 bent circulant: det", "rank %d" % qrank(C), "max hyperplane rank", max(r), "min", min(r))
# (3) n=4: every circulant has max hyperplane rank >= 4 (Turan argument); spot-check: all +-1 eigenvalue patterns
import numpy as np
best = 99
Hm = np.array([[(-1) ** bin(a & x).count("1") for x in range(16)] for a in range(16)])
for bits in range(1 << 15):           # d_0 = +1 fixed (scale), d_alpha = +-1
    d = np.array([1] + [(-1) ** (bits >> k & 1) for k in range(15)])
    mx = 0
    for h in range(1, 16):
        cnt = sum(1 for a in range(16) if a < (a ^ h) and d[a] + d[a ^ h] != 0)
        mx = max(mx, cnt)
    best = min(best, mx)
print("n=4 circulants with +-1 spectrum: min over spectra of max hyperplane rank =", best)
