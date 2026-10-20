#!/usr/bin/env python3
"""vgeom.py (v-fp-capacity): independent recomputation of fp-capacity's geometry (Lemma S, blocking sets,
rigid-pair rho and lambda). Own construction of the dual polar graph; nothing of fp-capacity or sl-* imported.
Symplectic form on F_2^(2m): <(a|b),(a'|b')> = a.b' + a'.b. Vector v = a | (b << m). L0 = {(0|b)}, L1 = {(b|b)}
(= Kbar L0, Kbar = prod of transvections T_(e_i|0): (a|b) -> (a + b | b)). Lagrangians found as the orbit of L0
under all symplectic transvections. usage: python3 vgeom.py m_exact_max samples_m4"""
import sys, itertools, random, json
import numpy as np
from collections import Counter

def run(m, sample=None, seed=7):
    N = 1 << (2*m); mask = (1 << m) - 1
    def form(u, v):
        a, b = u & mask, u >> m; a2, b2 = v & mask, v >> m
        return (bin(a & b2).count("1") + bin(a2 & b).count("1")) & 1
    def span(gens):
        S = {0}
        for g in gens: S |= {s ^ g for s in S}
        return frozenset(S)
    L0 = span([1 << (m + j) for j in range(m)])
    L1 = span([(1 << j) | (1 << (m + j)) for j in range(m)])
    # orbit under transvections
    trans = list(range(1, N))
    def T(v, L): return frozenset(x ^ (v if form(x, v) else 0) for x in L)
    seen = {L0}; frontier = [L0]
    while frontier:
        nf = []
        for L in frontier:
            for v in trans:
                L2 = T(v, L)
                if L2 not in seen: seen.add(L2); nf.append(L2)
        frontier = nf
    labs = sorted(seen, key=lambda s: sorted(s))
    V = len(labs)
    # sanity: each is totally isotropic of size 2^m
    assert all(len(L) == (1 << m) and all(form(x, y) == 0 for x in L for y in L) for L in labs[:50])
    B = np.zeros((V, N), dtype=np.int32)
    for i, L in enumerate(labs):
        for x in L: B[i, x] = 1
    inter = B @ B.T
    D = m - np.round(np.log2(inter)).astype(int)
    i0 = labs.index(L0); i1 = labs.index(L1)
    assert D[i0, i1] == m
    out = dict(m=m, frames=V)
    nb0 = np.flatnonzero(D[i0] == 1)
    out["neighbours_of_L0"] = len(nb0)
    # Lemma S: for F at level k, over neighbours N of L0, d(F,N)-k in {+1,0,-1}
    okS = True; prof = {}
    for F in range(V):
        k = int(D[i0, F]); d = D[F, nb0] - k
        a, z, dd = int((d == 1).sum()), int((d == 0).sum()), int((d == -1).sum())
        exp = (2*(2**m - 2**k), 2**k - 1, 2**k - 1)
        if (a, z, dd) != exp: okS = False
        prof.setdefault(k, Counter())[(a, z, dd)] += 1
        # characterisation: ascent iff U = F cap L0 not inside H = N cap L0
        U = labs[F] & L0
        for n_, Nn in zip(d, nb0):
            H = labs[Nn] & L0
            if (n_ == 1) != (not U <= H): okS = False
    out["lemmaS_ok"] = okS
    out["lemmaS_profile"] = {k: {str(kk): v for kk, v in c.items()} for k, c in prof.items()}
    # blocking sets: smallest set of frames at levels < m such that every neighbour N fails to ascend one of them
    if m <= 3 or sample:
        low = [F for F in range(V) if D[i0, F] < m]
        asc = {F: frozenset(int(n) for n in nb0 if D[F, n] == D[i0, F] + 1) for F in low}
        allnb = frozenset(int(n) for n in nb0)
        best = None
        for j in ((1, 2, 3) if m <= 3 else ()):
            for comb in itertools.combinations(low, j):
                inter_ = allnb
                for F in comb: inter_ = inter_ & asc[F]
                if not inter_: best = j; break
            if best: break
        out["min_blocking_set_levels_below_m"] = best
        # union-bound sufficiency: any bundle with sum(2^k_i - 1) < 2^m - 1 has a common ascent (check pairs+triples)
        bad = 0
        for j in ((2, 3) if m <= 3 else (2,)):
            for comb in itertools.combinations(low, j) if m <= 3 else itertools.islice(itertools.combinations(low, j), 200000):
                if sum(2**int(D[i0, F]) - 1 for F in comb) < 2**m - 1:
                    inter_ = allnb
                    for F in comb: inter_ = inter_ & asc[F]
                    if not inter_: bad += 1
        out["union_bound_violations"] = bad
    # rigid pairs
    f = D[i0] - D[i1]
    nopp = int((D[i0] == m).sum())
    out["ordered_opposite_pairs"] = V * nopp
    if sample:
        rng = random.Random(seed); pairs = []
        for _ in range(sample):
            A = rng.randrange(V); opps = np.flatnonzero(D[A] == m); pairs.append((A, int(opps[rng.randrange(len(opps))])))
    else:
        pairs = [(A, int(Ap)) for A in range(V) for Ap in np.flatnonzero(D[A] == m)]
    adj = (D == 1)
    rhoC = Counter(); lamC = Counter()
    for (A, Ap) in pairs:
        fA = D[A] - D[Ap]
        g2 = f + fA                    # = 2g, integer
        rho2 = int(g2.max() - g2.min())  # rho = rho2/4
        rhoC[rho2] += 1
        # longest path with each step f += 2 and fA += 2
        L = {}
        order = np.argsort(f, kind="stable")
        Lv = np.zeros(V, dtype=int)
        for X in order:
            for Y in np.flatnonzero(adj[X] & (f == f[X] + 2) & (fA == fA[X] + 2)):
                Lv[Y] = max(Lv[Y], Lv[X] + 1)
        lamC[int(Lv.max())] += 1
    n = len(pairs)
    out["pairs_examined"] = n
    out["rho_distribution"] = {str(k / 4): v for k, v in sorted(rhoC.items())}
    out["mean_rho"] = sum(k / 4 * v for k, v in rhoC.items()) / n
    out["count_rho_eq_m"] = rhoC.get(4 * m, 0)
    out["lambda_distribution"] = {str(k): v for k, v in sorted(lamC.items())}
    return out

if __name__ == "__main__":
    mmax = int(sys.argv[1]); s4 = int(sys.argv[2])
    res = []
    for m in range(2, mmax + 1):
        r = run(m); print(json.dumps(r), flush=True); res.append(r)
    if s4:
        r = run(4, sample=s4); r["sampled"] = True; print(json.dumps(r), flush=True); res.append(r)
    json.dump(res, open("vgeom_%s_%s.json" % (sys.argv[1], sys.argv[2]), "w"), indent=1)
