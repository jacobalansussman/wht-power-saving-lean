"""alsfam.py (d-R7): generate frame families (dims) passing octad Hall on P112 by annealing, then masked ALS
(read masks = Lines(F_k), serve free) with restarts; report best residual per family.
usage: alsfam.py dims nfam restarts iters seed [maxsecs]"""
import sys, time, json; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
from als import als
import os
dims = [int(x) for x in sys.argv[1].split(',')]; NF, NR, IT, SEED = map(int, sys.argv[2:6])
NCOORD = int(os.environ.get('NCOORD', '0'))     # fixed coordinate totals (dim 8, 56 ports) added in front
MAXS = float(sys.argv[6]) if len(sys.argv) > 6 else 100
t0 = time.time()
POOL = {d: pool(d, minports={6: 16, 7: 28, 8: 48}[d]) for d in set(dims)}
rng = np.random.default_rng(SEED)
COORD = np.array([[GRAM[o, P[T]] != 0 for T in range(NP)] for o in list(O0)])[:NCOORD]
def anneal(steps=20000):
    cur = [rng.integers(len(POOL[d][0])) for d in dims]
    S = np.vstack([COORD] + [POOL[d][1][i][None] for d, i in zip(dims, cur)])
    f = int(hall_deficit(S).clip(0).sum()); T = 3.0
    for step in range(steps):
        j = rng.integers(len(dims)); d = dims[j]; new = rng.integers(len(POOL[d][0]))
        S2 = S.copy(); S2[NCOORD + j] = POOL[d][1][new]
        f2 = int(hall_deficit(S2).clip(0).sum())
        if f2 <= f or rng.random() < np.exp((f - f2) / T): cur[j] = new; S = S2; f = f2
        T = max(0.05, T * 0.9995)
        if f == 0: break
    return f, cur, S
results = []
for fam in range(NF):
    if time.time() - t0 > MAXS: break
    f, cur, S = anneal()
    if f: print('family %d: Hall not reached (%d)' % (fam, f), flush=True); continue
    best = (1e9,)
    for rs in range(NR):
        res, X, R, h = als(S, ORTH, iters=IT, rng=rng)
        if res < best[0]: best = (res, X, R, h)
    res, X, R, h = best
    nt = natural_test(X @ R)
    duals_ = [POOL[d][0][i] for d, i in zip(dims, cur)]
    results.append({'dims': dims, 'duals': duals_, 'res': res})
    print('family %d: ports %s  best residual %.3e  (hist %s)  rank/nat %s  %.0fs' % (fam, S.sum(1).tolist(), res,
          ['%.1e' % x for x in h[::max(1, len(h) // 5)]], np.round(nt, 4), time.time() - t0), flush=True)
json.dump(results, open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/alsfam_c%d_%s_%d.json') % (NCOORD, '-'.join(map(str, dims)), SEED), 'w'))
print('BEST over families: %.3e' % min([r['res'] for r in results] + [1e9]))
