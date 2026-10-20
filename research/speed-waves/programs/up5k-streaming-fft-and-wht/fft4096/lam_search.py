# ap-6: search a uniform scale lam for the W16 rotator constants {1, c1, s1, c2} that minimises the total CSD digits
import itertools, math
import numpy as np
K = 14
vals = {}
for n in range(1, 5):
    for comb in itertools.combinations(range(0, K + 1), n):
        for sg in itertools.product((1, -1), repeat=n):
            if sg[0] != 1: continue
            v = sum(s * 2.0 ** -k for s, k in zip(sg, comb))
            if v <= 0: continue
            vals.setdefault(n, []).append((v, tuple(zip(sg, comb))))
for n in vals:
    vals[n].sort()
arr = {n: np.array([v for v, _ in vals[n]]) for n in vals}
def best(c, tol):
    for n in range(1, 5):
        a = arr[n]; i = np.searchsorted(a, c)
        for j in (i - 1, i):
            if 0 <= j < len(a) and abs(a[j] - c) <= tol:
                return n, vals[n][j]
    return None
C = [1.0, math.cos(math.pi / 8), math.sin(math.pi / 8), math.cos(math.pi / 4)]
res = []
for tol_db in (-72, -66):
    tol = 10 ** (tol_db / 20)
    bestl = None
    for lam in np.linspace(0.70, 1.0, 30001):
        tot = 0; reps = []
        ok = True
        for c in C:
            b = best(lam * c, tol * lam)
            if b is None: ok = False; break
            tot += b[0]; reps.append(b)
        if ok and (bestl is None or tot < bestl[0]):
            bestl = (tot, lam, reps)
    print('tol %d dB: total digits %d (baseline lam=1: 1+4+4+5=14), lam %.6f' % (tol_db, bestl[0], bestl[1]))
    for c, (n, (v, d)) in zip(('one', 'c1', 's1', 'c2'), bestl[2]):
        print('   %s: %d digits %s value %.7f (target %.7f)' % (c, n, d, v, bestl[1] * C[('one', 'c1', 's1', 'c2').index(c)]))
