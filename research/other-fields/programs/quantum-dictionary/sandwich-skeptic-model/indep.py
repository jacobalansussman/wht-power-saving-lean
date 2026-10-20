#!/usr/bin/env python3
"""indep.py (v-sandwich-model) -- independent checks, standard library only, no frames.py.
Vectors of F_2^(2m): (x, z) packed as v = x | z << m, omega = x.z' + x'.z.  L_0 = {(0|z)}, L_F = {(u|u)},
L_U = {(u | u+v): u in U, v in U^perp} (frames.py's convention, re-implemented).
1. all Lagrangians by brute force; interval lemma; excess histogram.
2. the m=3 witness: frames rebuilt from the printed bases; a symplectic phi found by my own backtracking;
   every membership and the timeline checked.
3. ordered pairs of distinct odd labels, orthogonal vs all, m = 2..9.
"""
import sys, itertools
def popc(x): return bin(x).count("1")
def span(vs):
    s = {0}
    for v in vs:
        if v not in s: s |= {x ^ v for x in s}
    return frozenset(s)
def om(u, v, m):
    M = (1 << m) - 1
    return (popc((u & M) & (v >> m)) + popc((v & M) & (u >> m))) & 1
def lagrangians(m):
    N = 2 * m; res = set()
    def rec(basis, S):
        if len(basis) == m:
            res.add(S); return
        for v in range(1, 1 << N):
            if v in S: continue
            if all(om(v, b, m) == 0 for b in basis) and v > max(basis, default=0):
                # canonical-ish: still may duplicate, set dedups
                rec(basis + [v], span(list(S) + [v]))
    rec([], frozenset([0]))
    return sorted(res, key=lambda s: sorted(s))
def perp(U, m): return frozenset(v for v in range(1 << m) if all(popc(u & v) % 2 == 0 for u in U))
def L_of(U, m): return frozenset(u | ((u ^ v) << m) for u in U for v in perp(U, m))
def dist(A, B, m): return m - (len(A & B).bit_length() - 1)

def part1(m):
    Ls = lagrangians(m); n = len(Ls)
    L0 = frozenset(z << m for z in range(1 << m)); LF = frozenset(z | (z << m) for z in range(1 << m))
    subs = {span(vs) for k in range(m + 1) for vs in itertools.combinations(range(1, 1 << m), k)}
    subframes = {L_of(U, m) for U in subs}
    I0F = {Z for Z in Ls if dist(L0, Z, m) + dist(Z, LF, m) == m}
    bad = 0; cnt = 0
    for X in Ls:
        for Y in Ls:
            d = dist(X, Y, m); cnt += 1
            for Z in Ls:
                inI = dist(X, Z, m) + dist(Z, Y, m) == d
                gen = span(list((Z & X) | (Z & Y))) == Z
                if inI != gen: bad += 1
    hist = {}
    for Z in Ls:
        e = dist(L0, Z, m) + dist(Z, LF, m) - m; hist[e] = hist.get(e, 0) + 1
    print("m=%d: %d Lagrangians; %d subspace frames; I(L0,LF)==subspace frames: %s; interval lemma on %d ordered pairs: %d failures; excess %s"
          % (m, n, len(subframes), I0F == subframes, cnt, bad, dict(sorted(hist.items()))))
    return Ls

def part2():
    m = 3
    def pv(s):   # printed "(xxx|zzz)", first char = bit 0
        xs, zs = s.strip("()").split("|")
        x = sum(1 << i for i, c in enumerate(xs) if c == "1"); z = sum(1 << i for i, c in enumerate(zs) if c == "1")
        return x | (z << m)
    def fr(*bs): return span([pv(b) for b in bs])
    F = {57: fr("(011|000)", "(100|100)", "(000|011)"), 134: fr("(111|100)", "(111|010)", "(111|001)"),
         71: fr("(000|100)", "(000|010)", "(000|001)"), 47: fr("(101|000)", "(000|010)", "(000|101)"),
         55: fr("(011|000)", "(000|100)", "(000|011)"), 51: fr("(101|000)", "(010|010)", "(000|101)"),
         79: fr("(100|100)", "(000|010)", "(000|001)")}
    for k, L in F.items():
        assert len(L) == 8 and all(om(a, b, m) == 0 for a in L for b in L), k
    L0 = frozenset(z << m for z in range(8)); LF = frozenset(z | (z << m) for z in range(8))
    assert F[71] == L0
    uS, uT = 0b001, 0b111
    LuS, LuT = L_of(span([uS]), m), L_of(span([uT]), m)
    LeS, LeT = L_of(perp(span([uS]), m), m), L_of(perp(span([uT]), m), m)
    d = lambda A, B: dist(A, B, m)
    inI = lambda Z, A, B: d(A, Z) + d(Z, B) == d(A, B)
    XS = lambda Z: inI(Z, LuS, LF); XT = lambda Z: inI(Z, LuT, LF)
    YS = lambda Z: inI(Z, L0, LeS); YT = lambda Z: inI(Z, L0, LeT)
    a0, ast, mu, nu = F[57], F[71], F[55], F[79]
    pa0, past, pmu, pnu = F[134], F[47], F[51], F[71]
    print("witness memberships: a0 in X_S %s, phi a0 in X_T %s, a* in Y_S %s, phi a* in Y_T %s, mu in Y_S %s, phi mu in X_T %s, nu in X_S %s, phi nu in Y_T %s"
          % (XS(a0), XT(pa0), YS(ast), YT(past), YS(mu), XT(pmu), XS(nu), YT(pnu)))
    print("distances: d(a0,a*)=%d d(a0,mu)=%d d(mu,a*)=%d | d(phi a0,phi mu)=%d d(phi mu,phi a*)=%d | old-style d(a0,phi a0)=%d d(a*,phi a*)=%d d(mu,phi mu)=%d d(nu,phi nu)=%d"
          % (d(a0, ast), d(a0, mu), d(mu, ast), d(pa0, pmu), d(pmu, past), d(a0, pa0), d(ast, past), d(mu, pmu), d(nu, pnu)))
    # timeline: X_S: nu(=<u_S>) -> a0 -> ... -> F ;  X_T: phi a0 (=<u_T>) -> phi mu -> F ;  R: a0 -> mu -> a*
    # Y_T: phi nu (=L0) -> phi a* -> ... -> e_T ;  Y_S: a* (=L0) -> mu -> ... -> e_S
    chain = lambda *Zs: sum(d(Zs[i], Zs[i + 1]) for i in range(len(Zs) - 1)) == d(Zs[0], Zs[-1])
    print("timeline geodesic: X_S <u_S>->nu->a0->F %s; X_T <u_T>->phi a0->phi mu->F %s; Y_S L0->a*->mu->e_S %s; Y_T L0->phi nu->phi a*->e_T %s; nu==<u_S> %s, phi a0==<u_T> %s"
          % (chain(LuS, nu, a0, LF), chain(LuT, pa0, pmu, LF), chain(L0, ast, mu, LeS), chain(L0, pnu, past, LeT), nu == LuS, pa0 == LuT))
    # find phi in Sp(6,2) with phi(a0)=pa0, phi(a*)=pa*, phi(mu)=pmu, phi(nu)=pnu.  f_i=(0|e_i) span L0 = a*, so f_i' in pa*.
    N = 6; e = [1 << i for i in range(m)]; f = [1 << (m + i) for i in range(m)]
    def apply(cols, v):
        r = 0
        for b in range(N):
            if (v >> b) & 1: r ^= cols[b]
        return r
    found = []
    cand_f = [v for v in past if v]
    def rec(k, imgs):
        if found: return
        if k == m:
            cols = [imgs[i][0] for i in range(m)] + [imgs[i][1] for i in range(m)]
            if all(frozenset(apply(cols, v) for v in A) == B for A, B in ((a0, pa0), (mu, pmu), (nu, pnu), (ast, past))):
                found.append(cols)
            return
        used = [v for p in imgs for v in p]
        for fi in cand_f:
            if any(om(fi, w, m) for w in used) or fi in span(used): continue
            for ei in range(1, 1 << N):
                if om(ei, fi, m) != 1 or any(om(ei, w, m) for w in used): continue
                if ei in span(used + [fi]): continue
                rec(k + 1, imgs + [(ei, fi)])
    rec(0, [])
    cols = found[0]
    ok = all(om(apply(cols, a), apply(cols, b), m) == om(a, b, m) for a in range(64) for b in range(64))
    print("phi found:", [format(c, "06b") for c in cols], "symplectic:", ok)
    # Gram floor check: p'_T, q'_S in R's coordinates (R's reference = x_S, so psi_S = id, psi_T = phi)
    inv = {apply(cols, v): v for v in range(64)}
    pS = uS | (uS << m); pT = uT | (uT << m); qS = uS << m; qT = uT << m   # vectors common to X_u frames / Y_u frames
    assert all(pS in Z for Z in [LuS, LF]) and all(qS in Z for Z in [L0, LeS])
    print("p_u = (u|u) lies in every X_u frame: %s; (u|0) in L_<u>: %s" % (pT in LuT and pT in LF, uT in LuT))
    pTp, qTp = inv[pT], inv[qT]
    print("R-coordinates: p'_S=p_S in a0 %s, p'_T in a0 %s, q'_S in a* %s, q'_T in a* %s; Gram [[w(p'_S,q'_S), w(p'_S,q'_T)],[w(p'_T,q'_S), w(p'_T,q'_T)]] = %s"
          % (pS in a0, pTp in a0, qS in ast, qTp in ast, [[om(pS, qS, m), om(pS, qTp, m)], [om(pTp, qS, m), om(pTp, qTp, m)]]))

def part3():
    for m in range(2, 10):
        odd = [u for u in range(1, 1 << m) if popc(u) & 1]
        allp = sum(1 for a in odd for b in odd if a != b)
        orth = sum(1 for a in odd for b in odd if a != b and popc(a & b) % 2 == 0)
        print("m=%d odd labels %d: ordered pairs %d (claim %d), orthogonal %d (claim %d), ratio %.4f"
              % (m, len(odd), allp, 2 ** (m - 1) * (2 ** (m - 1) - 1), orth, (2 ** (m - 1) - 1) * 2 ** (m - 2), allp / orth))

if __name__ == "__main__":
    what = sys.argv[1:] or ["1", "2", "3"]
    if "1" in what:
        part1(2)
        if "skip3" not in what: part1(3)
    if "2" in what: part2()
    if "3" in what: part3()
