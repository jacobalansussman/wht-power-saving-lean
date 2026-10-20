"""census.py m -- exact census of first backward steps from the kernel frame that lie on shortest routes, in the
widened game, by orbits of AGL(m,2) acting on perfect matchings of the 2^m rows of K (AGL(m,2) is the permutation
part of Lambda' = B cap K B K^-1, which maps routes to routes; checked on random images below)."""
import sys
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import sys, time, random, json
import numpy as np
from rt import Game, all_matchings, agl_perms

m = int(sys.argv[1])
g = Game(m); n = g.n
t0 = time.time()
G = agl_perms(m)
print("m=%d |AGL|=%d" % (m, len(G)), flush=True)
def gen_partners(n):
    tot = 1
    for k in range(n - 1, 0, -2):
        tot *= k
    P = np.zeros((tot, n), dtype=np.int8)
    cur = [-1] * n
    cnt = [0]
    def rec():
        try:
            p = cur.index(-1)
        except ValueError:
            P[cnt[0]] = cur; cnt[0] += 1; return
        for q in range(p + 1, n):
            if cur[q] == -1:
                cur[p] = q; cur[q] = p
                rec()
                cur[p] = -1; cur[q] = -1
    rec()
    assert cnt[0] == tot
    return P
partners = gen_partners(n)
print("matchings", len(partners), flush=True)
def to_matching(prow):
    return [(x, int(prow[x])) for x in range(n) if x < prow[x]]
class _L:
    def __len__(self): return len(partners)
    def __getitem__(self, i): return to_matching(partners[i])
allM = _L()
def enc_partner(P):  # P: (k, n) partner arrays -> uint64 codes
    P = P.astype(np.uint64)
    c = np.zeros(P.shape[0], dtype=np.uint64)
    for x in range(n):
        c |= P[:, x] << np.uint64(4 * x)
    return c

codes = enc_partner(partners)
order = np.argsort(codes); scodes = codes[order]
seen = np.zeros(len(allM), dtype=bool)
ar = np.arange(len(G))[:, None]

def orbit_codes(p):
    img = np.zeros_like(G)
    img[ar, G] = G[:, p.astype(np.int64)]          # img[i, G[i,x]] = G[i, p[x]]
    return np.unique(enc_partner(img))

E0 = g.compat(g.K)
assert len(E0) == n * (n - 1) // 2
orbits = []
for i in range(len(allM)):
    if seen[i]:
        continue
    oc = orbit_codes(partners[i])
    pos = order[np.searchsorted(scodes, oc)]
    assert (codes[pos] == oc).all()
    seen[pos] = True
    orbits.append((i, len(oc)))
print("orbits", len(orbits), "sum sizes", sum(s for _, s in orbits), "t=%.1f" % (time.time() - t0), flush=True)

res = []
rng = random.Random(1)
for (i, size) in orbits:
    M = allM[i]
    g3 = g.back_step(g.K, M, E0)
    ok = g.continues(g3, m - 1)
    nr = g.nroutes(g3, m - 1, {}) if ok else 0
    diffs = sorted(p ^ q for (p, q) in M)
    trans = len(set(diffs)) == 1
    # support structure of g3: bipartite row/col components
    S = np.abs(g3) > 1e-12
    supsets = sorted(set(tuple(np.nonzero(S[r])[0]) for r in range(n)))
    res.append(dict(rep=M, size=size, cont=ok, nroutes=nr, diffs=diffs, translation=trans,
                    n_distinct_supports=len(supsets)))
    print("orbit size %7d cont=%s routes=%d translation=%s diffs=%s nsupp=%d" % (size, ok, nr, trans, diffs, len(supsets)), flush=True)

# sanity: AGL images of continuing reps continue (random checks)
for r in res:
    M = r["rep"]
    for _ in range(3):
        s = G[rng.randrange(len(G))]
        M2 = [(int(s[p]), int(s[q])) for (p, q) in M]
        M2 = [(min(p, q), max(p, q)) for (p, q) in M2]
        g3 = g.back_step(g.K, M2, E0)
        assert g.continues(g3, m - 1) == r["cont"], "AGL action does not preserve continuation!"
print("AGL-invariance spot check passed")
tot = sum(r["size"] for r in res if r["cont"])
totr = sum(r["size"] * r["nroutes"] for r in res if r["cont"])
print("CONTINUING first steps: %d of %d; translation ones: %d; total routes (chains of cosets) = %d" % (
    tot, len(allM), sum(r["size"] for r in res if r["cont"] and r["translation"]), totr))
json.dump(res, open("census_m%d.json" % m, "w"), indent=1)
print("t=%.1f" % (time.time() - t0))
