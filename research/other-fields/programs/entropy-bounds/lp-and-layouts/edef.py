# edef.py -- cross deficiency e(D) of a three-pair layout: min over two of the pairs {i, j} and orientations (u of i,
# w of j; partners u', w') of (m - d(u, w')) + (m - d(w, u')). Lemma W+ (PROOF.md section 8): cut value >= 3m - e(D).
import sys, random, itertools, collections
import w3cut
def edef(D, m, P):
    best = 99
    for i, j in ((0, 1), (0, 2), (1, 2)):
        for (u, ub) in (P[i], P[i][::-1]):
            for (w, wb) in (P[j], P[j][::-1]):
                best = min(best, (m - int(D[u, wb])) + (m - int(D[w, ub])))
    return best
if __name__ == '__main__':
    m = int(sys.argv[1]); N = int(sys.argv[2]); seed = int(sys.argv[3])
    g = w3cut.G(m); rng = random.Random(seed)
    c = collections.Counter(edef(g.D, m, w3cut.rand_3k2(g, rng)) for _ in range(N))
    print(f"m={m}: e(D) distribution over {N} random layouts (six distinct frames): {dict(sorted(c.items()))}")
    # adversarial: hill-climb to maximise e(D) by moving terminals (keeping pairs opposite)
    worst = 0; wP = None
    for r in range(int(sys.argv[4]) if len(sys.argv) > 4 else 0):
        P = w3cut.rand_3k2(g, rng); e = edef(g.D, m, P)
        for it in range(400):
            k = rng.randrange(3); end = rng.randrange(2)
            s, t = P[k]; x = (s, t)[end]; y = (t, s)[end]
            nb = rng.choice(g.adj[x])
            if g.D[nb, y] != m or nb in {v for p in P for v in p}: continue
            Q = list(P); Q[k] = (nb, y) if end == 0 else (y, nb)
            e2 = edef(g.D, m, Q)
            if e2 >= e: P, e = Q, e2
        if e > worst: worst, wP = e, P
    if wP: print(f"m={m}: hill-climb max e(D) = {worst} at {wP}")
