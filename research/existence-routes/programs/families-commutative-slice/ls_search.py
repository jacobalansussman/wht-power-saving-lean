"""ls_search.py: local search for commutative-slice move multisets that PASS the deletion (two-point) test with
M < m W, then the tight-pair clique test (Theorem 1) on each one found.
usage: ls_search.py m W M tries seed"""
import sys, time
import numpy as np
from refl import quad_forms
from clique_test import clique_test

m, W, M, tries, seed = map(int, sys.argv[1:6])
vals, tr = quad_forms(m)
V = vals[tr == 1][:, 1:].astype(np.int32)        # odd forms x directions (1..2^m-1)
nd = V.shape[1]
rng = np.random.default_rng(seed)


def score(w):
    mn = w.min()
    return (min(mn, W), -int((w == mn).sum()) if mn < W else 0)


found = {}
t0 = time.time()
for t in range(tries):
    c = np.zeros(nd, dtype=np.int32)
    for d in rng.integers(0, nd, size=M):
        c[d] += 1
    w = V @ c
    best = score(w)
    for step in range(int(sys.argv[6]) if len(sys.argv) > 6 else 600):
        if best[0] >= W:
            break
        # try all single transfers a -> b; take the best (random tie-break)
        srcs = np.nonzero(c)[0]
        bestc = None
        for a in srcs:
            Wab = (w - V[:, a])[:, None] + V            # forms x b
            mins = Wab.min(axis=0)
            cnts = (Wab == mins[None, :]).sum(axis=0)
            for b in range(nd):
                if b == a:
                    continue
                s = (min(int(mins[b]), W), -int(cnts[b]) if mins[b] < W else 0)
                key = (s, rng.random())
                if bestc is None or key > bestc[0]:
                    bestc = (key, a, b)
        (s, _), a, b = bestc
        if s <= best and rng.random() < 0.3:
            a = rng.choice(srcs); b = int(rng.integers(0, nd))
            s = score(w - V[:, a] + V[:, b])
        c[a] -= 1; c[b] += 1; w = w - V[:, a] + V[:, b]; best = s
    if best[0] >= W:
        key = tuple(c.tolist())
        if key not in found:
            S = [d + 1 for d in range(nd) for _ in range(c[d])]
            found[key] = S
            v = clique_test(m, S, W, verbose=False)
            print(f"try {t}: PASSES deletion: S={S} -> clique test {v}", flush=True)
print(f"m={m} W={W} M={M}: tries {tries}, distinct deletion-passers {len(found)}, t={time.time()-t0:.0f}s", flush=True)
