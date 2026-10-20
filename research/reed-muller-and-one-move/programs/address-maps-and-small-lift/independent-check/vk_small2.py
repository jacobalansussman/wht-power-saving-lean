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

m, W, M = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ci, nc = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (0, 1)
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
    """max weighted clique; adj[v] = bitmask of neighbours (over direction indices)."""
    cands = list(cands); best = 0
    def rec(Pmask, s):
        nonlocal best
        if Pmask == 0:
            if s > best: best = s
            return
        tot = 0; mm = Pmask
        while mm:
            b = mm & -mm; tot += w[b.bit_length() - 1]; mm ^= b
        if s + tot <= best: return
        while Pmask:
            b = Pmask & -Pmask; v = b.bit_length() - 1
            rec(Pmask & adj[v], s + w[v])
            Pmask ^= b
            tot -= w[v]
            if s + tot <= best: return
    m0 = 0
    for v in cands: m0 |= 1 << v
    rec(m0, 0)
    return best

def tests_CZ(cnt, tl):
    """cnt: counts per direction (len q-1); tl: list of target indices. returns 'C', 'Z' or None (survives)."""
    D = (Qtab.astype(np.int16) @ cnt.astype(np.int16))
    r = P[tl].sum(0).astype(np.int16)
    tight = np.nonzero((D == r) & (r >= 1))[0]
    if len(tight) == 0: return None
    Qt = Qtab[tight].astype(np.int32) * (cnt > 0)[None, :]
    livev = Qt.any(0)
    if not livev.any(): return None
    J = (Qt.T @ Qt) > 0
    adj = [0] * (q - 1)
    for a in range(q - 1):
        if livev[a]:
            row = 0
            for b in np.nonzero(J[a] & livev)[0]:
                if b != a: row |= 1 << int(b)
            adj[a] = row
    wv = [int(c) for c in cnt]
    live = [int(z) for z in np.nonzero(livev)[0]]
    c = max_clique_w(live, adj, wv)
    if c > W: return 'C'
    need = np.nonzero(D + r < 2 * c)[0]
    for F in need:
        sub = [z for z in live if Qtab[F, z]]
        if not sub: continue
        y = max_clique_w(sub, adj, wv)
        if D[F] + r[F] < 2 * y: return 'Z'
    return None

t0 = time.time()
stats = dict(multisets=0, A_ok=0, del_pass=0, C=0, Z=0, surv=0)
surv = []
idx = -1
for comb in itertools.combinations_with_replacement(range(q - 1), M):
    idx += 1
    if idx % nc != ci: continue
    stats['multisets'] += 1
    cnt = np.bincount(np.array(comb), minlength=q - 1)
    D = (Qtab.astype(np.int16) @ cnt.astype(np.int16))
    A = D - pI
    if (A < 0).any(): continue
    stats['A_ok'] += 1
    zero = (A == 0)
    cand = np.where(~(P[:, zero].any(1)))[0]           # T with p_T <= A
    if len(cand) == 0: continue
    # recursive search over multisets of W-1 targets from cand (nondecreasing)
    Pc = P[cand].astype(np.int16)
    def rec(start, chosen, Arem):
        if len(chosen) == W - 1:
            stats['del_pass'] += 1
            tl = [ti] + chosen
            v = tests_CZ(cnt, tl)
            if v is None:
                stats['surv'] += 1; surv.append((comb, tuple(tl)))
            else:
                stats[v] += 1
            return
        ok = np.nonzero(~((Pc[start:] > Arem[None, :]).any(1)))[0] + start
        for k in ok:
            rec(int(k), chosen + [int(cand[k])], Arem - Pc[k])
    rec(0, [], A)
print(f"m={m} W={W} M={M} chunk {ci}/{nc}: targets {nT}; " + ", ".join(f"{k} {v}" for k, v in stats.items()) +
      f"; {time.time()-t0:.0f}s", flush=True)
for s in surv[:40]:
    print("SURV", s[0], [tuple(targets[t].flatten()) for t in s[1]])
