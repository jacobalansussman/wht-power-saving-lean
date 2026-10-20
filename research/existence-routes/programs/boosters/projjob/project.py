# b56-methods project.py: the PROJECTION RECURSION.  For a point p of F_2^{2m} (symplectic) and a Lagrangian L,
# pi_p(L) = (L cap p^perp + <p>)/<p>, a Lagrangian of p^perp/p (rank m-1).  Given any family Y' at rank m-1,
# Y_L = (+)_p Y'_p(pi_p L) is a family at rank m and abar_r = sum_s N_r(s) abar'_s, where N_r(s) = number of
# points p with dist(pi_p L, pi_p M) = s for one (any) pair L, M at distance r (Sp acts distance-transitively,
# and pi_p is Stab(p)-equivariant, so the image of the pairs at distance r is uniform within each distance s).
# This file counts N_r(s) exactly (all points, one pair per r), and checks the count on random pairs.
# Usage: python3 -I project.py m
import sys, random
m = int(sys.argv[1]); n2 = 2 * m
# vectors: (x|z) as int with x in low m bits, z in high m bits; form = x.z' + x'.z
def form(a, b):
    ax, az = a & ((1 << m) - 1), a >> m; bx, bz = b & ((1 << m) - 1), b >> m
    return (bin(ax & bz).count('1') + bin(bx & az).count('1')) & 1
def span(B):
    out = [0]
    for b in B:
        if b in out: continue
        out = out + [x ^ b for x in out]
    return out
def dim(B): return len(span(B)).bit_length() - 1
def lag(r):  # Z-frame on first m-r qubits, X-frame on last r qubits; L0 = Z-frame
    B = []
    for q in range(m):
        if q < m - r: B.append(1 << (m + q))   # z_q
        else: B.append(1 << q)                 # x_q
    return B
L0 = span(lag(0))
def proj(Lspan, p):
    # (L cap p^perp + p): return as set of vectors (contains p); the quotient intersection dims are computed upstairs
    S = [v for v in Lspan if form(v, p) == 0]
    return set(span(S + [p]))
def dist_proj(Ls, Ms, p):
    A = proj(Ls, p); B = proj(Ms, p)
    inter = len(A & B).bit_length() - 1   # dim of intersection upstairs (contains p)
    return (m - 1) - (inter - 1)
from collections import Counter
res = {}
for r in range(1, m + 1):
    M = span(lag(r)); assert len(set(L0) & set(M)) == 2 ** (m - r)
    c = Counter(dist_proj(L0, M, p) for p in range(1, 1 << n2))
    res[r] = dict(sorted(c.items()))
    print(f"m={m} r={r}: N_r(s) = {res[r]}  (total points {sum(c.values())})")
# random check: apply a random symplectic map (product of random transvections) to the pair, counts must agree
rng = random.Random(1)
def transv(v, t):  # T_t(v) = v + <v,t> t
    return v ^ t if form(v, t) else v
for trial in range(3):
    ts = [rng.randrange(1, 1 << n2) for _ in range(6 * m)]
    def g(v):
        for t in ts: v = transv(v, t)
        return v
    for r in (1, m):
        L = span([g(v) for v in lag(0)]); M = span([g(v) for v in lag(r)])
        c = Counter(dist_proj(L, M, p) for p in range(1, 1 << n2))
        assert dict(sorted(c.items())) == res[r], (r, c)
print("random-pair checks agree")
