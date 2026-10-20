# w3compat.py -- the l-infinity-plane construction for three opposite pairs.
# Pairs i, j are COMPATIBLE if for some choice of ends u of pair i and w of pair j (partners u', w'):
#     d(u, w') = m and d(w, u') = m.
# Then phi1 = min(d(u, .), d(w, .)) is 1-Lipschitz and stretches pairs i and j to m, phi2 = d(x, .) stretches the
# third pair (x one of its ends), and v -> (phi1(v), phi2(v)) in the l-infinity plane (isometric to the l1 plane)
# is a 1-Lipschitz l1 contraction exact on all three pairs: cut value 3m.
import sys, random, numpy as np, itertools
import w3cut
def compatible(D, m, P, i, j):
    for (u, ub) in (P[i], P[i][::-1]):
        for (w, wb) in (P[j], P[j][::-1]):
            if D[u, wb] == m and D[w, ub] == m: return (u, w)
    return None
def any_compat(D, m, P):
    for i, j in ((0, 1), (0, 2), (1, 2)):
        c = compatible(D, m, P, i, j)
        if c: return (i, j, c)
    return None
if __name__ == '__main__':
    m = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3])
    g = w3cut.G(m); rng = random.Random(seed); D = g.D
    cnt = 0; bad = []
    for tr in range(trials):
        P = w3cut.rand_3k2(g, rng)
        if any_compat(D, m, P): cnt += 1
        else: bad.append(P)
    print(f"m={m}: {cnt}/{trials} layouts have a compatible pair of pairs")
    for P in bad[:4]:
        X = [v for p in P for v in p]
        print("  none:", P, [[int(D[a, b]) for b in X] for a in X])
    import json; json.dump(bad[:200], open(f"w3compat_bad_m{m}_s{seed}.json", "w"))
