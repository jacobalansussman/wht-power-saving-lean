"""v-rm2-phase: own class-based D + C + Z enumeration at block 4, any W.
usage: v_cls.py W M mode [shard nshards]
  mode 'alt': target W-multisets from all 448 with at least one alternating; 'noalt': only the 420 congruent to I.
Direction multisets up to GL(4,2) (own orbit computation); targets: every nondecreasing W-tuple (alternating indexed
last, so 'at least one alternating' = the last index is alternating), deletion-pruned depth first; then cliques and
the zero-block test at the leaves (rm-theory 8.5, 8.7)."""
import sys, itertools, time, numpy as np
from vb import *

m = 4; W, M, mode = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
sh, nsh = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (0, 1)
t0 = time.time()
Qz = Qz_matrix(m); nf, nd = Qz.shape
Tall = sym_nonsing(m)
na = [t for t in Tall if t.diagonal().any()]; al = [t for t in Tall if not t.diagonal().any()]
T = na + (al if mode == 'alt' else []); NA = len(na); nt = len(T)
PT = pair_matrix(m, T)                       # nf x nt
G = gl(m)
perm = np.zeros((len(G), nd), dtype=np.int64)
for g, A in enumerate(G):
    for z in range(1, 2 ** m):
        w = (A @ vec(z, m)) % 2
        perm[g, z - 1] = int(sum(int(w[i]) << i for i in range(m))) - 1
pw = (16 ** np.arange(M)).astype(np.int64)
seen = set(); classes = []
for ms in itertools.combinations_with_replacement(range(nd), M):
    k = int(np.dot(ms, pw)) if False else sum(z * 16 ** i for i, z in enumerate(ms))
    if k in seen: continue
    im = np.sort(perm[:, list(ms)], axis=1) @ pw
    seen.update(im.tolist()); classes.append(ms)
del seen
print(f"[{W} in {M} mode={mode} shard {sh}/{nsh}] {len(classes)} direction classes ({time.time()-t0:.0f} s)", flush=True)

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

st = dict(leaves=0, C=0, Z=0, surv=0); survs = []
for ci, ms in enumerate(classes):
    if ci % nsh != sh: continue
    cnt = np.zeros(nd, dtype=np.int64)
    for d in ms: cnt[d] += 1
    Dq = cnt @ Qz.T
    inst = list(ms); n = len(inst)
    fm = [sum(1 << k for k, d in enumerate(inst) if Qz[f, d]) for f in range(nf)]
    def leaf(r):
        tight = [fm[f] for f in np.where((Dq == r) & (r >= 1))[0]]
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
        for f in np.where(Dq + r < 2 * c)[0]:
            if Dq[f] + r[f] < 2 * mc(fm[f] & live): return 'Z'
        return None
    def dfs(tup, s, depth):
        lo = tup[-1] if tup else 0
        cand = np.arange(lo, nt)
        ok = cand[(PT[:, lo:] <= s[:, None]).all(0)]
        if depth == W - 1 and mode == 'alt': ok = ok[ok >= NA]
        for j in ok:
            s2 = s - PT[:, j]
            if depth == W - 1:
                st['leaves'] += 1
                v = leaf(Dq - s2)
                if v is None:
                    st['surv'] += 1
                    if len(survs) < 50: survs.append((ms, tup + (int(j),)))
                else: st[v] += 1
            else:
                dfs(tup + (int(j),), s2, depth + 1)
    dfs((), Dq.copy(), 0)
print(f"[{W} in {M} mode={mode} shard {sh}/{nsh}] done: leaves passing deletion {st['leaves']}, killed by C {st['C']}, "
      f"killed by Z {st['Z']}, survivors {st['surv']} ({time.time()-t0:.0f} s)")
for x in survs[:20]: print("  survivor", x)
