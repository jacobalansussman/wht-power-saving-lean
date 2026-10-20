"""zpool.py <zb dump> <k list>: pool the zero blocks E[Y,Y] = 0 derived by the star recursion (global facts about E),
then (a) clique test: a set Y with E[Y,Y] = 0 and |Y| > W is impossible (K[Y,Y] = I + Gamma_Y unit triangular,
rank K <= W: ex2-families Theorem 1's argument); (b) matching test: for every form F, rank E[D_F,D_F] = s_F needs a
matching of size s_F among the cells of D_F x D_F not forced to zero (Corollary 10')."""
import sys, struct, itertools
import numpy as np
zb = [l.split() for l in open(sys.argv[1]) if l.startswith("ZB")]
ks = [int(x) for x in sys.argv[2].split(",")]; W = len(ks)
Z = set()
for p in zb:
    Y = [z for z in range(1, 32) if int(p[4]) >> (z - 1) & 1]
    for x in Y:
        for y in Y: Z.add((x, y))
print("forced zero cells", len(Z), "zero diagonal", sum(1 for x in range(1, 32) if (x, x) in Z))
# (a) clique: vertices with zero diagonal, edges x~y if both (x,y),(y,x) zero
V = [x for x in range(1, 32) if (x, x) in Z]
adj = {x: {y for y in V if y != x and (x, y) in Z and (y, x) in Z} for x in V}
best = []
def bk(R, P, X):
    global best
    if not P and not X:
        if len(R) > len(best): best = list(R)
        return
    for v in list(P):
        bk(R + [v], P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}
bk([], set(V), set())
print("largest zero clique", len(best), best, "W", W, "VIOLATION" if len(best) > W else "")
# (b) matchings
raw = open("tables5.bin", "rb").read(); nf = struct.unpack_from("<i", raw, 0)[0]; o = 4
mask = np.frombuffer(raw, "<u4", nf, o).astype(np.int64); o += 4 * nf
wt = np.frombuffer(raw, "<i4", nf, o).astype(np.int64); o += 4 * nf
vid = np.frombuffer(raw, "<i4", nf, o).astype(np.int64); o += 4 * nf; o += 4
pk = np.frombuffer(raw, "<u4", 32, o).astype(np.int64)
cnt = np.zeros(31, dtype=np.int64)
for k in ks: cnt[k % 31] += 1
w = np.array([sum(cnt[k] for k in range(31) if pk[v] >> k & 1) for v in range(32)]); s = wt - w[vid]
def maxmatch(D):
    match = {}
    def aug(x, seen):
        for y in D:
            if (x, y) in Z or y in seen: continue
            seen.add(y)
            if y not in match or aug(match[y], seen): match[y] = x; return True
        return False
    return sum(1 for x in D if aug(x, set()))
bad = 0; tight = 0
for F in range(1, nf):
    D = [z for z in range(1, 32) if int(mask[F]) >> (z - 1) & 1]
    if not any((x, y) in Z for x in D for y in D): continue
    mm = maxmatch(D)
    if mm < s[F]: bad += 1; print("MATCHING VIOLATION F", F, "matching", mm, "< s_F", s[F]) if bad < 5 else None
    elif mm == s[F]: tight += 1
print("matching violations", bad, "forms at equality", tight)
