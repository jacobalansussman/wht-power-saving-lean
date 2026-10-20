# Refinement for the near-uniform cases nearU6.py left open (e3-slice6): star recursion with MIXED bounds.
# For an even 4-flat F avoiding Z: sigma_A = min slack of odd 4-flats G (avoiding 0) through each hyperplane A of F;
# planes: (g-2) k <= sum sigma - 2|Psi| with g = 3; lines: g = 7 planes; points: g = 15 lines.  Entries E_xy can be
# nonzero only if the line {x,y} (or, for x = y, the point) allows rank >= 1; if the allowed pattern has no perfect
# matching, E[F,F] is singular, contradicting Theorem 10(b).
import itertools, sys
exec(open("nearU6.py").read().split("dirs = range(1, N)")[0])
def flats_in(Fl, d):
    pts = sorted(Fl); a = pts[0]; D = [p ^ a for p in pts]
    subs = set()
    for vs in itertools.combinations([x for x in D if x], d):
        S = frozenset(span(vs))
        if len(S) == 1 << d: subs.add(S)
    out = set()
    for S in subs:
        for t in D:
            out.add(frozenset(p for p, dd in zip(pts, D) if (dd ^ t) in S))
    return list(out)
def has_perfect_matching(nodes, allowed):
    match = {}
    def try_(u, seen):
        for v in nodes:
            if allowed(u, v) and v not in seen:
                seen.add(v)
                if v not in match or try_(match[v], seen):
                    match[v] = u; return True
        return False
    return all(try_(u, set()) for u in nodes)
def refined(Z, W):
    Zs = list(Z)
    for Fl in even4:
        if any(z in Fl for z in Zs): continue
        sig = {}
        ok = True
        for A in H_of[Fl]:
            cands = [16 + sum(1 for z in Zs if z in G) - W for G in odd_through.get(A, []) if True]
            if A not in odd_through: cands = [G for G in odd4 if A <= G]; cands = [16 + sum(1 for z in Zs if z in G) - W for G in cands]
            sig[A] = min(cands + [8])
        planes = flats_in(Fl, 2); lines = flats_in(Fl, 1)
        sp = {}
        for P in planes:
            Hs = [A for A in H_of[Fl] if P <= A]
            sp[P] = min(4, sum(sig[A] for A in Hs) - 8)          # g = 3, |Psi| = 4
        sl = {}
        for L in lines:
            Ps = [P for P in planes if L <= P]
            sl[L] = min(2, (sum(sp[P] for P in Ps) - 4) // 5)    # g = 7, |Psi| = 2
        if min(sp.values()) < 0 or min(sl.values()) < 0: return True
        spt = {}
        for x in Fl:
            Ls = [L for L in lines if x in L]
            spt[x] = min(1, (sum(sl[L] for L in Ls) - 2) // 13)  # g = 15, |Psi| = 1
        if min(spt.values()) < 0: return True
        nodes = sorted(Fl)
        lineof = {frozenset(L): L for L in lines}
        allowed = lambda u, v: (spt[u] >= 1) if u == v else (sl[frozenset([u, v])] >= 1)
        if not has_perfect_matching(nodes, allowed): return True
    return False
cases = [((63,), 13)] + [(Z, 13) for Z in itertools.combinations_with_replacement(range(1, 64), 2)] + [((63, 63), 14)]
left = []
for Z, W in cases:
    if excluded(Z, W): continue
    if refined(Z, W): print("now excluded:", Z, "W =", W, flush=True)
    else: left.append((Z, W))
print("still open:", len(left), left[:10])
