"""contfam.py (d-R7): masked decoders by continuation ALS (penalty on off-mask entries rising to a hard mask).
Families: 'coord' (8 coordinate frames, control), or saved Hall-passing families from out/*.json, or fresh annealing.
usage: contfam.py coord|dims nfam restarts iters seed maxsecs"""
import sys, time, json, os; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
from als import als_cont
kind = sys.argv[1]; NF, NR, IT, SEED = map(int, sys.argv[2:6]); MAXS = float(sys.argv[6])
rng = np.random.default_rng(SEED); t0 = time.time()
COORD = np.array([[GRAM[o, P[T]] != 0 for T in range(NP)] for o in list(O0)])
def fams():
    if kind == 'coord':
        yield 'coord8', COORD
        yield 'coord7+one7', None
        return
    dims = [int(x) for x in kind.split(',')]
    POOL = {d: pool(d, minports={6: 16, 7: 28, 8: 48}[d]) for d in set(dims)}
    for fam in range(NF):
        cur = [rng.integers(len(POOL[d][0])) for d in dims]
        S = np.array([POOL[d][1][i] for d, i in zip(dims, cur)]); f = int(hall_deficit(S).clip(0).sum()); T = 3.0
        for step in range(20000):
            j = rng.integers(len(dims)); d = dims[j]; new = rng.integers(len(POOL[d][0]))
            S2 = S.copy(); S2[j] = POOL[d][1][new]; f2 = int(hall_deficit(S2).clip(0).sum())
            if f2 <= f or rng.random() < np.exp((f - f2) / T): cur[j] = new; S = S2; f = f2
            T = max(0.05, T * 0.9995)
            if f == 0: break
        if f == 0: yield 'fam%d' % fam, S
        else: print('fam%d: Hall not reached (%d)' % (fam, f), flush=True)
for name, S in fams():
    if time.time() - t0 > MAXS: break
    if S is None: continue
    best = None
    for rs in range(NR):
        res, X, R, h = als_cont(S, ORTH, iters=IT, rng=rng)
        if best is None or res < best[0]: best = (res, X, R, h)
    res, X, R, h = best
    print('%s: ports %s  best residual after continuation %.3e  path %s  nat %s  %.0fs' % (name, S.sum(1).tolist(), res,
          [(a, '%.1e' % b) for a, b in h[::max(1, len(h) // 6)]] + [h[-1]], np.round(natural_test(X @ (R * S)), 4), time.time() - t0), flush=True)
