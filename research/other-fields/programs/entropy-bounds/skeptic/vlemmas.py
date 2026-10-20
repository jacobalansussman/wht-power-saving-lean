# vlemmas.py -- v-nc-gain: direct checks of nc-gain session 2's lemmas on my own Gamma_m.
#  (1) index map onto nc-gain's dpg.py (subspace equality), edge sets and distances equal;
#  (2) near-polygon claim (b) of Lemma W+, brute force (all triples at m = 2, 3; sampled opposite pairs at m = 4);
#  (3) Lemma W: for random two-class stretchable layouts, the explicit packing
#        delta(u,v) = |f1(u)-f1(v)| + |f2(u)-f2(v)|, f1 = (phiA+phiB)/2, f2 = (phiA-phiB)/2
#      has delta <= 1 on every edge and delta(s_j,t_j) = m on every pair (exact, in halves);
#  (4) Lemma W+: for random / hill-climbed three-pair layouts, build D' by the moves of claim (b), check D' compatible,
#      and check the Lemma W packing of D' gives sum_j delta'(s_j,t_j) >= 3m - e(D) on D itself;
#  (5) e(D) distribution on random layouts (six distinct frames) at m = 2, 3, 4.
import sys, random, itertools, collections, time
import numpy as np
from fractions import Fraction
sys.path.insert(0, '.')
import vgamma

def phi_of(g, lows):
    return g.D[lows].min(axis=0).astype(np.int64)

def stretch_orient(g, P, idx):
    m = g.m
    for o in itertools.product((0, 1), repeat=len(idx)):
        lo = [P[i][o[k]] for k, i in enumerate(idx)]; hi = [P[i][1 - o[k]] for k, i in enumerate(idx)]
        if all(g.D[a, b] == m for a in lo for b in hi): return lo, hi   # includes a==b pairs (opposite)
    return None

def two_class(g, P):
    W = len(P)
    for mask in range(1 << W):
        A = [i for i in range(W) if not mask >> i & 1]; B = [i for i in range(W) if mask >> i & 1]
        if not A: continue
        a = stretch_orient(g, P, A)
        if a is None: continue
        b = stretch_orient(g, P, B) if B else ([], [])
        if b is not None: return A, a, B, b
    return None

def lemmaW_packing(g, P, cls):
    A, (loA, hiA), B, (loB, hiB) = cls
    pA = phi_of(g, loA)
    pB = phi_of(g, loB) if loB else np.zeros(g.n, np.int64)
    # work in doubled units: 2 f1 = pA + pB, 2 f2 = pA - pB
    F1 = pA + pB; F2 = pA - pB
    eu = np.array([u for u, v in g.edges]); ev = np.array([v for u, v in g.edges])
    load2 = np.abs(F1[eu] - F1[ev]) + np.abs(F2[eu] - F2[ev])      # = 2 * delta(e)
    pairs2 = [abs(int(F1[s] - F1[t])) + abs(int(F2[s] - F2[t])) for s, t in P]
    return int(load2.max()), pairs2, (F1, F2)

def edef(g, P):
    m = g.m; best = 99; arg = None
    for i, j in itertools.combinations(range(len(P)), 2):
        for (u, up) in (P[i], P[i][::-1]):
            for (w, wp) in (P[j], P[j][::-1]):
                val = (m - int(g.D[u, wp])) + (m - int(g.D[w, up]))
                if val < best: best, arg = val, (i, j, u, up, w, wp)
    return best, arg

def step_b(g, u, y, w):
    """claim (b): y opposite u, d(w,y) = k < m -> a neighbour y' of y opposite u with d(w,y') = k+1."""
    m = g.m
    for z in g.adj[y]:
        if g.D[u, z] == m and g.D[w, z] == g.D[w, y] + 1: return int(z)
    return None

def rand_layout(g, rng, W):
    m = g.m
    while True:
        P = []
        for _ in range(W):
            s = rng.randrange(g.n); opp = np.nonzero(g.D[s] == m)[0]; P.append((s, int(rng.choice(list(opp)))))
        if len({x for p in P for x in p}) == 2 * W: return P

def main(m, quick=False):
    t0 = time.time(); rng = random.Random(1000 + m)
    g = vgamma.VG(m); n = g.n
    out = []
    # (1) map to dpg
    sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'lp-and-layouts'))
    import dpg
    verts, idx, adj, dist = dpg.build(m) if m < 4 else dpg_build_fast(m)
    perm = np.array([idx[S] for S in g.L])           # my index -> dpg index
    assert sorted(perm.tolist()) == list(range(n))
    myE = {frozenset((int(perm[u]), int(perm[v]))) for u, v in g.edges}
    thE = {frozenset((u, v)) for u in range(n) for v in adj[u]}
    Dth = np.array(dist, np.int16) if dist is not None else None
    sameD = bool((g.D == Dth[np.ix_(perm, perm)]).all()) if Dth is not None else 'skipped'
    out.append(f"(1) m={m}: dpg map ok; edge sets equal: {myE == thE}; distances equal: {sameD}")
    print(out[-1], flush=True)
    # (2) near-polygon claim (b)
    bad = 0; tested = 0
    if m <= 3:
        pairs = [(u, y) for u in range(n) for y in np.nonzero(g.D[u] == m)[0]]
    else:
        pairs = []
        while len(pairs) < 300:
            u = rng.randrange(n); y = int(rng.choice(list(np.nonzero(g.D[u] == m)[0]))); pairs.append((u, y))
    for (u, y) in pairs:
        y = int(y)
        nb = g.adj[y]; oppu = nb[g.D[u, nb] == m]
        for w_ in range(n):
            k = g.D[w_, y]
            if k >= m: continue
            tested += 1
            if not (g.D[w_, oppu] == k + 1).any(): bad += 1
    out.append(f"(2) m={m}: claim (b) tested on {tested} triples (u opposite y, d(w,y) < m; {len(pairs)} (u,y) pairs): failures {bad}")
    print(out[-1], flush=True)
    # (3) Lemma W explicit packing
    for W in (3, 4, 5, 6):
        ntot = 0; nst = 0; nbad = 0; T = 300 if not quick else 60
        for _ in range(T):
            P = rand_layout(g, rng, W); ntot += 1
            cls = two_class(g, P)
            if cls is None: continue
            nst += 1
            mx, pairs2, _ = lemmaW_packing(g, P, cls)
            if mx > 2 or any(p != 2 * m for p in pairs2): nbad += 1
        out.append(f"(3) m={m} W={W}: {nst}/{ntot} random layouts two-class stretchable; explicit Lemma W packing fails on {nbad}")
        print(out[-1], flush=True)
    # (4) Lemma W+ explicit construction
    def wplus(P):
        e, (i, j, u, up, w, wp) = edef(g, P)
        a = m - int(g.D[w, up]); b = m - int(g.D[u, wp])
        y = up; steps1 = 0
        while g.D[w, y] < m:
            y2 = step_b(g, u, y, w)
            if y2 is None: return e, 'stuck1', None
            y = y2; steps1 += 1
        z = wp; steps2 = 0
        while g.D[u, z] < m:
            z2 = step_b(g, w, z, u)
            if z2 is None: return e, 'stuck2', None
            z = z2; steps2 += 1
        assert g.D[w, y] == m and g.D[u, z] == m and g.D[u, y] == m and g.D[w, z] == m
        Q = list(P); Q[i] = (u, y); Q[j] = (w, z)
        cls = two_class(g, Q)
        if cls is None: return e, 'Dprime not two-class', None
        mx, pairs2, (F1, F2) = lemmaW_packing(g, Q, cls)
        valD = sum(abs(int(F1[s] - F1[t])) + abs(int(F2[s] - F2[t])) for s, t in P)   # doubled
        okload = mx <= 2
        return e, (steps1 == a and steps2 == b, okload, Fraction(valD, 2)), Q
    cnt = collections.Counter(); fails = 0; checked = 0
    for _ in range(3000 if not quick else 500):
        P = rand_layout(g, rng, 3)
        e, r, Q = wplus(P); cnt[e] += 1
        if e == 0: continue
        checked += 1
        if not isinstance(r, tuple) or not r[0] or not r[1] or r[2] < 3 * m - e: fails += 1; print('  W+ FAIL', P, e, r)
    out.append(f"(4,5) m={m}: e(D) over {sum(cnt.values())} random layouts {dict(sorted(cnt.items()))} "
               f"(shares {', '.join(f'{k}: {100*v/sum(cnt.values()):.2f}%' for k, v in sorted(cnt.items()))}); "
               f"Lemma W+ construction checked on {checked} with e >= 1: failures {fails}")
    print(out[-1], flush=True)
    # hill-climb to large e and check W+ there too
    worst = 0; hc = 0; hfail = 0
    for r_ in range(40 if not quick else 10):
        P = rand_layout(g, rng, 3); e = edef(g, P)[0]
        for it in range(400):
            k = rng.randrange(3); end = rng.randrange(2); s, t = P[k]; x = (s, t)[end]; y = (t, s)[end]
            nb = int(rng.choice(list(g.adj[x])))
            if g.D[nb, y] != m or nb in {v for p in P for v in p}: continue
            Q = list(P); Q[k] = (nb, y) if end == 0 else (y, nb)
            e2 = edef(g, Q)[0]
            if e2 >= e: P, e = Q, e2
        worst = max(worst, e)
        e, r, _ = wplus(P); hc += 1
        if e > 0 and (not isinstance(r, tuple) or not r[0] or not r[1] or r[2] < 3 * m - e): hfail += 1; print('  W+ FAIL (hill)', P, e, r)
    out.append(f"(4b) m={m}: hill-climbed {hc} layouts, max e(D) = {worst}; Lemma W+ construction failures {hfail}  ({time.time()-t0:.0f}s)")
    print(out[-1], flush=True)
    return out

def dpg_build_fast(m):
    """dpg.build(m) but with the O(n^2) Python distance loop replaced (only verts/idx/adj are taken from dpg)."""
    import dpg, collections as C
    L0 = dpg.span([1 << (m + i) for i in range(m)])
    verts = [L0]; idx = {L0: 0}; adj = [[]]; q = C.deque([0])
    while q:
        i = q.popleft(); L = verts[i]; basis = dpg.basis_of(L); nb = set()
        for H in dpg.hyperplanes(L, basis, m):
            Hb = dpg.basis_of(H)
            perp = [v for v in range(1 << (2 * m)) if all(dpg.form(v, h, m) == 0 for h in Hb)]
            for v in perp:
                if v in H: continue
                N = frozenset(H | {x ^ v for x in H})
                if N != L: nb.add(N)
        for N in nb:
            if N not in idx: idx[N] = len(verts); verts.append(N); adj.append([]); q.append(idx[N])
            adj[i].append(idx[N])
    for i in range(len(verts)): adj[i] = sorted(set(adj[i]))
    return verts, idx, adj, None

if __name__ == '__main__':
    m = int(sys.argv[1]); quick = len(sys.argv) > 2
    main(m, quick)
