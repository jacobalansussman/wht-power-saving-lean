"""v-rm-theory: independent small-lift enumerator (rm-theory section 5b), own code and own enumeration order.
Instead of direction classes up to GL with all target tuples, it FIXES T_1 = I (GL is transitive on targets
congruent to I) and enumerates EVERY direction multiset (no GL reduction, no lemma filter) and every multiset of
the other W-1 targets. Tests, all exact (from rm-theory 8.5, 8.7, which I re-derived):
  (D) deletion incl. span: for every form Q, sum_l <Q,T_l> <= |D_Q| (|D_Q| = moves z with Q(z) = 1)
  (C) cliques: live directions (in some tight support), joined if in a common tight support; max weighted clique <= W
  (Z) zero block: for every form F and every clique Y inside D_F: |D_F| + r_F >= 2|Y|
usage: vk_small.py m W M [chunk_index nchunks]
"""
import itertools, sys, time
import numpy as np

m = int(sys.argv[2]); W = 4
q = 1 << m
pts = list(range(1, q))
pairs = [(i, j) for i in range(m) for j in range(i, m)]
nf = 1 << len(pairs)
Qbits = np.array([[(Q >> t) & 1 for t in range(len(pairs))] for Q in range(nf)], dtype=np.int64)
X = np.array([[(z >> j) & 1 for j in range(m)] for z in pts], dtype=np.int64)       # (q-1) x m
mon = np.array([X[:, i] * X[:, j] for (i, j) in pairs]).T                             # (q-1) x npairs
Qtab = ((Qbits @ mon.T) % 2).astype(np.int8)                                         # nf x (q-1): Q(z)
def pvec(T):
    tv = np.array([T[i, j] for (i, j) in pairs], dtype=np.int64)
    return ((Qbits @ tv) % 2).astype(np.int8)
# targets congruent to I = g g^T, g invertible: generate from GL directly (own route, not via det/diag test)
def gl_elems():
    cols = pts
    for c in itertools.product(cols, repeat=m):
        G = np.array([[(c[j] >> i) & 1 for j in range(m)] for i in range(m)], dtype=np.int64)
        # invertible over F2?
        A = G.copy() % 2; r = 0
        for col in range(m):
            pv = [i for i in range(r, m) if A[i, col]]
            if not pv: continue
            A[[r, pv[0]]] = A[[pv[0], r]]
            for i in range(m):
                if i != r and A[i, col]: A[i] ^= A[r]
            r += 1
        if r == m: yield G
tset = {}
for G in gl_elems():
    T = (G @ G.T) % 2
    key = tuple(T.flatten())
    if key not in tset: tset[key] = T
targets = list(tset.values())
I_key = tuple(np.eye(m, dtype=np.int64).flatten())
ti = [k for k in tset].index(I_key)
P = np.array([pvec(T) for T in targets])                                              # nT x nf
pI = P[ti].astype(np.int16)
nT = len(targets)

def max_clique_w(cands, adj, w):
    best = 0
    def rec(R, Pset, s):
        nonlocal best
        if s + sum(w[v] for v in Pset) <= best: return
        if not Pset:
            best = max(best, s); return
        for v in list(Pset):
            rec(R + [v], [u for u in Pset if u != v and adj[v][u]], s + w[v])
            Pset = [u for u in Pset if u != v]
            if s + sum(w[x] for x in Pset) <= best: return
    rec([], list(cands), 0)
    return best

def tests_CZ(cnt, tl):
    """cnt: counts per direction (len q-1); tl: list of target indices. returns 'C', 'Z' or None (survives)."""
    D = (Qtab.astype(np.int16) @ cnt.astype(np.int16))
    r = P[tl].sum(0).astype(np.int16)
    tight = np.where((D == r) & (r >= 1))[0]
    live = [z for z in range(q - 1) if cnt[z] > 0 and Qtab[tight, z].any()] if len(tight) else []
    if not live: return None
    adj = {a: {b: bool((Qtab[tight, a] & Qtab[tight, b]).any()) for b in live} for a in live}
    wv = {z: int(cnt[z]) for z in live}
    c = max_clique_w(live, adj, wv)
    if c > W: return 'C'
    # zero block
    need = np.where(D + r < 2 * c)[0]
    for F in need:
        sub = [z for z in live if Qtab[F, z]]
        if not sub: continue
        y = max_clique_w(sub, adj, wv)
        if D[F] + r[F] < 2 * y: return 'Z'
    return None


# targets in rm-theory's t_small order (bits over pairs ascending, nonzero diagonal, invertible): map to my index
def rank2(T):
    A = T.copy() % 2; r = 0
    for c in range(m):
        pv = [i for i in range(r, m) if A[i, c]]
        if not pv: continue
        A[[r, pv[0]]] = A[[pv[0], r]]
        for i in range(m):
            if i != r and A[i, c]: A[i] ^= A[r]
        r += 1
    return r
theirs = []
for bits in range(1 << len(pairs)):
    T = np.zeros((m, m), dtype=np.int64)
    for t, (i, j) in enumerate(pairs):
        if (bits >> t) & 1: T[i, j] = T[j, i] = 1
    if not T.diagonal().any() or rank2(T) != m: continue
    theirs.append(tuple(T.flatten()))
keys = list(tset.keys())
assert sorted(theirs) == sorted(keys), (len(theirs), len(keys))
mymap = [keys.index(k) for k in theirs]
fn = sys.argv[1]
res = {}
for line in open(fn):
    line = line.strip()
    if not line or '|' not in line: continue
    d, t = line.split('|')
    dirs = [int(x) for x in d.split(',')]; tl = [mymap[int(x)] for x in t.split(',')]
    W = len(tl)
    cnt = np.bincount(np.array(dirs) - 1, minlength=q - 1)
    D = (Qtab.astype(np.int16) @ cnt.astype(np.int16)); r = P[tl].sum(0)
    if (r > D).any(): v = 'D'
    else: v = tests_CZ(cnt, tl) or 'survives'
    res[v] = res.get(v, 0) + 1
print(fn, res)
