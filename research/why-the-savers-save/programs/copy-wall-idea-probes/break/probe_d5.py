# SUPERSEDED by probe_petersen.py (pure-Python gradient too slow; killed, no output)
"""probe_d5.py (g-break, H2) -- D5 bit unit: is there a fitting decoder of rank 4 (= alpha) over R / Q, beating the
natural F_2 rank 5?  Fitting: B[S,T] = 0 when lines S,T share exactly one coordinate (non-orthogonal), B[S,S] != 0.
Try B = Gram of vectors x_S in R^4 with x_S = x_{support(S)} (both signs share a vector; same-support lines know
each other).  Then need x_A . x_B = 0 whenever supports A, B in C(5,2) share one index: an orthonormal
representation of the Petersen graph in R^4.  Also: exact F_2 minrank of the D5 pattern by DFS (is it 5?)."""
import itertools, random, math, sys
pairs = list(itertools.combinations(range(5), 2))
cons = [(a, b) for a in range(10) for b in range(a + 1, 10) if len(set(pairs[a]) & set(pairs[b])) == 1]
def try_dim(d, seeds=200):
    best = 1e9
    for sd in range(seeds):
        rng = random.Random(sd)
        X = [[rng.gauss(0, 1) for _ in range(d)] for _ in range(10)]
        for it in range(4000):
            for x in X:
                n = math.sqrt(sum(t * t for t in x)); x[:] = [t / n for t in x]
            # gradient step on sum (x_a.x_b)^2
            grad = [[0.0] * d for _ in range(10)]
            err = 0.0
            for a, b in cons:
                s = sum(p * q for p, q in zip(X[a], X[b])); err += s * s
                for k in range(d):
                    grad[a][k] += 2 * s * X[b][k]; grad[b][k] += 2 * s * X[a][k]
            if err < 1e-20: break
            for a in range(10):
                X[a] = [x - 0.1 * g for x, g in zip(X[a], grad[a])]
        best = min(best, err)
        if err < 1e-16:
            return err, X
    return best, None
for d in (3, 4, 5):
    err, X = try_dim(d, seeds=40)
    print('Petersen orthonormal representation in R^%d: residual %.3g %s' % (d, err, '(found)' if X else ''))
    if X:
        for p, x in zip(pairs, X):
            print('   ', p, ['%+.4f' % t for t in x])
        break
# exact F_2 minrank of the D5 pattern
lines = []
for (i, j) in pairs:
    for s in (1, -1):
        v = [0] * 5; v[i] = 1; v[j] = s; lines.append(v)
n = len(lines)
Gm = [[sum(a * b for a, b in zip(x, y)) for y in lines] for x in lines]
forb = [[T for T in range(n) if T != S and Gm[S][T] != 0] for S in range(n)]
known = [[T for T in range(n) if T != S and Gm[S][T] == 0] for S in range(n)]
def reduce(x, basis):
    for b in basis:
        x = min(x, x ^ b)
    return x
def minrank_f2(r):
    order = list(range(n))
    def dfs(k, basis):
        if k == n: return True
        S = order[k]
        base = 1 << S
        kn = known[S]
        for m in range(1 << len(kn)):
            x = base
            for t, T in enumerate(kn):
                if m >> t & 1: x |= 1 << T
            y = reduce(x, basis)
            if y == 0:
                if dfs(k + 1, basis): return True
            elif len(basis) < r:
                if dfs(k + 1, sorted(basis + [y], reverse=True)): return True
        return False
    return dfs(0, [])
for r in (4, 5):
    print('D5 pattern: F_2 minrank <= %d ? %s' % (r, minrank_f2(r)))
