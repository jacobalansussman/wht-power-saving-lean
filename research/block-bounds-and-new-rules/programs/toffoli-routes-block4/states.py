"""states.py m -- the route states of the widened game: S_j = frames at depth j on shortest routes from B to BK.
S_{m-1} = continuing first backward steps (all of them, expanded from census orbit reps by AGL); S_{m-2} = their
continuing halving children. Counts and block structure (row/column components) of S_{m-2}. By the transpose map
Bh -> B h^-T K, |S_j| = |S_{m-j}|."""
import sys
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import json, time, collections, hashlib
import numpy as np
from rt import Game, agl_perms

m = int(sys.argv[1])
g = Game(m); n = g.n
cen = json.load(open("census_m%d.json" % m)  )
G = agl_perms(m)
E0 = g.compat(g.K)

def ncomp(h):
    S = np.abs(h) > 1e-12
    parent = list(range(2 * n))
    def f(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for r in range(n):
        for c in np.nonzero(S[r])[0]:
            a, b = f(r), f(n + int(c))
            if a != b:
                parent[a] = b
    roots = collections.Counter(f(x) for x in range(2 * n))
    return tuple(sorted(v // 2 for v in roots.values()))

t0 = time.time()
S3 = 0
S2 = {}
for r in cen:
    if not r["cont"]:
        continue
    M = r["rep"]
    seen = set()
    for s in G:
        M2 = tuple(sorted(tuple(sorted((int(s[p]), int(s[q])))) for (p, q) in M))
        if M2 in seen:
            continue
        seen.add(M2)
        h3 = g.back_step(g.K, list(M2), E0)
        for kk, (MM, c) in g.children(h3).items():
            if g.continues(c, m - 2):
                hk = hashlib.blake2b(kk, digest_size=16).digest()
                if hk not in S2:
                    S2[hk] = (ncomp(c), r["translation"])
    assert len(seen) == r["size"]
    S3 += len(seen)
    g.memo.clear()
    print("orbit size %d done; |S2| so far %d; t=%.1f" % (r["size"], len(S2), time.time() - t0), flush=True)
print("|S_%d| = %d, |S_%d| = %d" % (m - 1, S3, m - 2, len(S2)))
cnt = collections.Counter(v[0] for v in S2.values())
for k, v in sorted(cnt.items()):
    print("S_%d block sizes (row/col component sizes) %s : %d" % (m - 2, k, v))
