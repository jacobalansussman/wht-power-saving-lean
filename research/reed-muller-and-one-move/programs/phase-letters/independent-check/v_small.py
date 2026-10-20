"""v-rm2-phase: own exhaustive D + C + Z enumeration at block 4 (rm-theory 8.5-8.7 definitions).
usage: v_small.py W M first pool [chunk_index nchunks]
  first: 'alt' (T_1 = a fixed alternating target A0) or 'I' (T_1 = I)
  pool:  'all' (448 nonsingular symmetric) or 'noalt' (the 420 congruent to I)
Every direction multiset of size M (no GL reduction, no free-direction filter); T_2 <= ... <= T_W from the pool.
Prints: multisets, tuples passing deletion, killed by cliques, killed by zero-block, survivors (and lists survivors).
"""
import sys, itertools, numpy as np
from vb import *

m = 4
W, M = int(sys.argv[1]), int(sys.argv[2]); first, pool = sys.argv[3], sys.argv[4]
ci, nc = (int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else (0, 1)
Qz = Qz_matrix(m)                      # 1024 x 15
T = sym_nonsing(m)
if pool == 'noalt': T = [t for t in T if t.diagonal().any()]
alt = [t for t in sym_nonsing(m) if not t.diagonal().any()]
T1 = alt[0] if first == 'alt' else np.eye(m, dtype=np.int64)
PT = pair_matrix(m, T)                 # 1024 x nt
p1 = pair_matrix(m, [T1])[:, 0]
nf, nd = Qz.shape

def maxclique(cand, adj):
    best = 0
    def bk(R, P, X):
        nonlocal best
        if P == 0 and X == 0:
            best = max(best, R); return
        if R + bin(P).count('1') <= best: return
        u = (P | X) & -(P | X); ui = u.bit_length() - 1
        Q = P & ~adj[ui]
        while Q:
            v = Q & -Q; vi = v.bit_length() - 1
            bk(R + 1, P & adj[vi], X & adj[vi]); P &= ~v; X |= v; Q &= ~v
    bk(0, cand, 0)
    return best

def cz(cnt, tl):
    """cnt: count vector over 15 directions; tl: list of target indices (None = T1). Returns 'C', 'Z' or None."""
    r = p1.copy()
    for t in tl: r = r + PT[:, t]
    inst = [d for d in range(nd) for _ in range(cnt[d])]
    n = len(inst)
    dmask = np.zeros(nf, dtype=object)
    fm = []
    for f in range(nf):
        mk = 0
        for k, d in enumerate(inst):
            if Qz[f, d]: mk |= 1 << k
        fm.append(mk)
    Dq = cnt @ Qz.T
    tight = [fm[f] for f in range(nf) if Dq[f] == r[f] and r[f] >= 1]
    adj = [0] * n; live = 0
    for mk in tight:
        live |= mk
        for k in range(n):
            if mk >> k & 1: adj[k] |= mk & ~(1 << k)
    cache = {}
    def mc(mask):
        if mask not in cache: cache[mask] = maxclique(mask, adj)
        return cache[mask]
    c = mc(live)
    if c > W: return 'C'
    for f in range(nf):
        if Dq[f] + r[f] < 2 * c:
            if Dq[f] + r[f] < 2 * mc(fm[f] & live): return 'Z'
    return None

ms = list(itertools.combinations_with_replacement(range(nd), M))
ms = ms[ci::nc]
nms = len(ms); passD = []; CH = 50000
for c0 in range(0, nms, CH):
    blk = ms[c0:c0 + CH]
    cnt = np.zeros((len(blk), nd), dtype=np.int64)
    for i, t in enumerate(blk):
        for d in t: cnt[i, d] += 1
    Dq = cnt @ Qz.T
    s = Dq - p1[None, :]
    ok = np.where(s.min(1) >= 0)[0]
    for i in ok:
        si = s[i]
        if W == 1:
            passD.append((cnt[i], [])); continue
        allowed = np.where(PT[si == 0].sum(0) == 0)[0]
        if W == 2:
            for t in allowed: passD.append((cnt[i], [t]))
            continue
        if W == 3:
            A = PT[si == 1][:, allowed]          # forms with slack 1
            conf = A.T @ A                        # number of slack-1 forms where both odd
            ii, jj = np.where(np.triu(conf == 0))
            for a, b in zip(ii, jj): passD.append((cnt[i], [allowed[a], allowed[b]]))
            continue
        raise SystemExit("W > 3 not implemented")
kC = kZ = 0; surv = []
for cnt, tl in passD:
    v = cz(cnt, tl)
    if v == 'C': kC += 1
    elif v == 'Z': kZ += 1
    else: surv.append((tuple(int(x) for x in cnt), [int(t) for t in tl]))
print(f"block 4, {W} in {M}, T1={first}, pool={pool} ({len(T)} targets), chunk {ci}/{nc}: multisets {nms}, "
      f"pass deletion {len(passD)}, killed by C {kC}, killed by Z {kZ}, survivors {len(surv)}")
for x in surv[:20]: print("  survivor", x)
