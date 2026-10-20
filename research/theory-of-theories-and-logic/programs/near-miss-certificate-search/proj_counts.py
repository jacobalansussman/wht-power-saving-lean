# proj_counts.py -- cert-07 session 3: projection (Harish-Chandra) families from block 3 to block 4.
# For a point p of PG(7,2) and a frame (Lagrangian) U of F_2^8, pi_p(U) = ((U cap p^perp) + p)/p is a frame of the
# block-3 residue p^perp/p.  Given a block-3 family Z (profile psi3_d'), Y_U = (+)_p Z^(p)_{pi_p U} is a block-4 family
# with psi4_d = sum_p psi3_{d_p}, d_p = dist(pi_p U, pi_p V).  This script counts N(d, d') = #{p : d_p = d'} [exact].
import itertools
m = 4; n2 = 2 * m
def form(x, y):  # standard symplectic form on F_2^8: pairs (i, i+4)
    s = 0
    for i in range(m): s ^= ((x >> i) & 1) & ((y >> (i + m)) & 1) ^ ((x >> (i + m)) & 1) & ((y >> i) & 1)
    return s
def span(vs):
    S = {0}
    for v in vs: S |= {s ^ v for s in S}
    return frozenset(S)
U = span([1 << (i + m) for i in range(m)])          # Z-type Lagrangian
# a frame at distance d from U: keep the first m-d Z directions, replace the others by X directions
def frame_at(d): return span([1 << (i + m) for i in range(m - d)] + [1 << i for i in range(m - d, m)])
def dim(S): return len(S).bit_length() - 1
def proj_dist(p, A, B):
    Ap = span([v for v in A if form(v, p) == 0] + [p]); Bp = span([v for v in B if form(v, p) == 0] + [p])
    return (m - 1) - (dim(Ap & Bp) - 1)
for d in range(1, m + 1):
    V = frame_at(d); cnt = {}
    for p in range(1, 2 ** n2):
        dp = proj_dist(p, U, V); cnt[dp] = cnt.get(dp, 0) + 1
    print("d =", d, "N(d, d') =", dict(sorted(cnt.items())))
