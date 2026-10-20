# interval.py -- the "interval relaxation" Q_m of the best Theorem-R bound (pf-duality, 13 Oct 2026).
#
# The frames between two opposite frames L, M at block size m are N_s = s + (s^perp cap M), one for each
# subspace s of L = F_2^m; their distance is d(s,t) = dim(s+t) - dim(s cap t) (PROOF.md, Lemma I1).
# A FAMILY gives each s a subspace Y_s of one vector space; delta(s,t) = dim(Y_s+Y_t) - dim Y_t.
# ratio(Y) = [mean delta over ordered pairs at distance m] / [mean delta over ordered pairs at distance 1].
# Q_m = sup of ratio over all families over a field of characteristic 0 (PROOF.md, Theorem I).
# This file evaluates ratio for explicit families (exact rational ranks):
#   gf2      : Y_s = N_s over GF(2) (perfect: ratio m)
#   pointlift: Y_s = Q-span of fixed vectors v_x, x in s \ 0  (v_x = the 0/1 vector of x, or given)
#   cut      : Y_s = Q if s in S else 0   (coordinate families; exhaustive at m = 3)
import sys, itertools
from fractions import Fraction

def subspaces(m):
    """all subspaces of F_2^m as frozensets of vectors (ints), including {0}"""
    n = 1 << m
    seen = set(); out = []
    def close(gens):
        S = {0}
        for g in gens:
            S |= {x ^ g for x in S}
        return frozenset(S)
    # build by adding generators
    frontier = [frozenset({0})]
    seen.add(frozenset({0}))
    while frontier:
        new = []
        for S in frontier:
            for v in range(1, n):
                if v not in S:
                    T = frozenset(S | {x ^ v for x in S})
                    if T not in seen:
                        seen.add(T); new.append(T)
        frontier = new
    out = sorted(seen, key=lambda S: (len(S), sorted(S)))
    return out

def dimof(S):
    return len(S).bit_length() - 1

def structure(m):
    subs = subspaces(m)
    idx = {S: i for i, S in enumerate(subs)}
    dims = [dimof(S) for S in subs]
    N = len(subs)
    dist = [[0] * N for _ in range(N)]
    for i, S in enumerate(subs):
        for j, T in enumerate(subs):
            inter = S & T
            # dim(S+T) = dim S + dim T - dim(S cap T)
            dsum = dims[i] + dims[j] - dimof(inter)
            dist[i][j] = dsum - dimof(inter)
    return subs, dims, dist

def rank_mod(rows, p=2305843009213693951):
    """rank of a list of integer vectors modulo a large prime (= rational rank unless p is unlucky)"""
    rows = [list(r) for r in rows if any(r)]
    if not rows: return 0
    k = len(rows[0]); r = 0
    M = [[x % p for x in row] for row in rows]
    for c in range(k):
        piv = None
        for i in range(r, len(M)):
            if M[i][c]: piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
        if r == len(M): break
    return r

def rank_frac(rows):
    rows = [[Fraction(x) for x in row] for row in rows if any(row)]
    if not rows: return 0
    k = len(rows[0]); r = 0
    for c in range(k):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0: piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r

def evaluate(m, delta, subs, dist):
    """mean delta per distance class over ordered pairs; returns list abar[0..m]"""
    N = len(subs)
    tot = [Fraction(0)] * (m + 1); cnt = [0] * (m + 1)
    for i in range(N):
        for j in range(N):
            d = dist[i][j]
            tot[d] += delta(i, j); cnt[d] += 1
    return [tot[d] / cnt[d] for d in range(m + 1)], cnt

def pointlift_delta(subs, vec, exact=False):
    rk = rank_frac if exact else rank_mod
    span = [[vec[x] for x in sorted(S) if x] for S in subs]
    r1 = [rk(s) for s in span]
    cache = {}
    def delta(i, j):
        key = (i, j)
        if key not in cache:
            cache[key] = rk(span[i] + span[j]) - r1[j]
        return cache[key]
    return delta, r1

def bits(x, m):
    return [(x >> b) & 1 for b in range(m)]

if __name__ == "__main__":
    for m in [2, 3, 4]:
        subs, dims, dist = structure(m)
        N = len(subs)
        print(f"== m = {m}: {N} frames in the interval; dims histogram", {d: dims.count(d) for d in range(m + 1)})
        gf2 = lambda i, j: dist[i][j]
        ab, cnt = evaluate(m, gf2, subs, dist)
        print("  ordered pairs per distance:", cnt)
        print("  gf2 family: abar =", [str(a) for a in ab], " ratio abar_m/abar_1 =", ab[m] / ab[1])
        vec = {x: bits(x, m) for x in range(1 << m)}
        delta, r1 = pointlift_delta(subs, vec, exact=(m <= 3))
        ab, _ = evaluate(m, delta, subs, dist)
        print("  0/1 point-lift: abar =", [f"{float(a):.4f}" for a in ab], " ratio =", ab[m] / ab[1], "=", float(ab[m] / ab[1]),
              " abar_2/abar_1 =", float(ab[2] / ab[1]) if m >= 2 else None)
        sys.stdout.flush()
