# w3ident.py -- test the identification criterion: for three opposite pairs, is there a choice of two terminals u, w
# of different pairs such that identifying u ~ w keeps all three pair distances = m? Then the demand graph becomes
# P3 + K2 (inside two stars), the cut condition suffices there (Papernov-Lomonosov), and pulling back gives cut
# value 3m on Gamma_m.
import sys, random, numpy as np, itertools
import w3cut
m = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3])
g = w3cut.G(m); rng = random.Random(seed); D = g.D
def ok_ident(P):
    hits = []
    for j, k in itertools.permutations(range(3), 2):
        if j > k: continue
        x, xb = P[3 - j - k]
        for (u, ub) in (P[j], P[j][::-1]):
            for (w, wb) in (P[k], P[k][::-1]):
                if D[w, ub] == m and D[u, wb] == m and D[x, u] + D[w, xb] >= m and D[x, w] + D[u, xb] >= m:
                    hits.append((u, w))
    return hits
cnt = 0; bad = []
for tr in range(trials):
    P = w3cut.rand_3k2(g, rng)
    h = ok_ident(P)
    if h: cnt += 1
    else: bad.append(P)
print(f"m={m}: {cnt}/{trials} layouts have a good identification")
for P in bad[:5]:
    X = [v for p in P for v in p]
    print("  no identification:", P, "terminal distances:", [[int(D[a, b]) for b in X] for a in X])
