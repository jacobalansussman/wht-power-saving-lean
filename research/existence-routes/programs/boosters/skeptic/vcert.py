# v-e3-boost vcert.py: my own certificate profiles (A, C_k, distance sets on the dual polar graph), validated by brute
# force at m = 2, 3; my own walk condition; the survivors S6, S7, S8 and their stated properties.
from fractions import Fraction as Fr
from itertools import combinations, product
import math, sys

# ---------- dual polar graph of Sp(2m,2): intersection array (standard: b_i = q^(i+1)[m-i], c_i = [i], q = 2) ----------
def iarray(m):
    b = [2 ** (i + 1) * (2 ** (m - i) - 1) for i in range(m + 1)]
    c = [2 ** i - 1 for i in range(m + 1)]
    k = b[0]; a = [k - b[i] - c[i] for i in range(m + 1)]
    return a, b, c, k
def pmat(m):
    """p[d][i][j] = #{W : d(U,W) = i, d(V,W) = j} for d(U,V) = d."""
    a, b, c, k = iarray(m); n = m + 1
    L = [[0] * n for _ in range(n)]          # A1 * A_j = sum_i L[i][j] A_i
    for j in range(n):
        if j - 1 >= 0: L[j - 1][j] = b[j - 1]
        L[j][j] = a[j]
        if j + 1 < n: L[j + 1][j] = c[j + 1]
    def mul(X, Y): return [[sum(X[i][t] * Y[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    I = [[int(i == j) for j in range(n)] for i in range(n)]
    V = [I, L]                                # v_i(L), A_{i+1} = ((A1 - a_i) A_i - b_{i-1} A_{i-1}) / c_{i+1}
    for i in range(1, m):
        X = mul(L, V[i]); Y = [[X[r][s] - a[i] * V[i][r][s] - b[i - 1] * V[i - 1][r][s] for s in range(n)] for r in range(n)]
        assert all(y % c[i + 1] == 0 for row in Y for y in row)
        V.append([[y // c[i + 1] for y in row] for row in Y])
    # p^d_{ij} = (v_i(L))[d][j]
    p = [[[V[i][d][j] for j in range(n)] for i in range(n)] for d in range(n)]
    ks = [p[0][i][i] for i in range(n)]
    return p, ks
def dist_profile(m, S):
    p, ks = pmat(m)
    tot = sum(ks[i] for i in S)
    return [tot - sum(p[d][i][j] for i in S for j in S) for d in range(1, m + 1)]

# ---------- brute force at m = 2, 3: Lagrangians of F_2^{2m}, distances, and the distance-set profiles ----------
def brute(m):
    def omega(u, v):
        ux, uz, vx, vz = u >> m, u & ((1 << m) - 1), v >> m, v & ((1 << m) - 1)
        return (bin(ux & vz).count("1") + bin(uz & vx).count("1")) & 1
    vecs = range(1, 1 << (2 * m))
    lags = set()
    def span(gens):
        s = {0}
        for g in gens: s |= {x ^ g for x in s}
        return frozenset(s)
    # enumerate isotropic bases greedily (all m-subsets would be too many at m=3? 63 choose 3 = 39711: fine)
    for gens in combinations(vecs, m):
        if all(omega(a, b) == 0 for a, b in combinations(gens, 2)):
            s = span(gens)
            if len(s) == 2 ** m: lags.add(s)
    lags = list(lags)
    L0 = span([1 << i for i in range(m)])        # z-part = low bits: the diagonal (Z-type) frame
    def dist(U, V): return m - int(math.log2(len(U & V)))
    prof = {}
    dists = [[dist(U, V) for V in lags] for U in lags]
    idx0 = lags.index(L0)
    for r in range(1, m + 2):
        for S in combinations(range(m + 1), r):
            SU = lambda u: {w for w in range(len(lags)) if dists[u][w] in S}
            S0 = SU(idx0)
            pr = []
            for d in range(1, m + 1):
                V = next(v for v in range(len(lags)) if dists[idx0][v] == d)
                pr.append(len(SU(V) - S0))
            prof[S] = pr
    return len(lags), prof

if __name__ == "__main__":
    for m in (2, 3):
        n, prof = brute(m)
        ok = all(prof[S] == dist_profile(m, S) for S in prof)
        a, b, c, k = iarray(m)
        print(f"m={m}: {n} frames (expect {math.prod(2**i + 1 for i in range(1, m+1))}); all {len(prof)} distance-set profiles"
              f" match the intersection-number formula: {ok}; S={{0,1}}: {dist_profile(m,(0,1))}")
    for m in range(2, 9):
        p, ks = pmat(m)
        print(f"m={m}: sphere sizes {ks}, total {sum(ks)} = frames {math.prod(2**i + 1 for i in range(1, m+1))}")
