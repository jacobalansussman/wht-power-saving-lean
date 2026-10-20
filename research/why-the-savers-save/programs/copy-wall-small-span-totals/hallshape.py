"""hallshape.py (d-R7): annealing for frame multisets (any dims 1..9) that pass octad Hall on P112.
Reports the least total Hall deficit found per shape.  usage: hallshape.py 'dims;dims;...' seeds steps maxsecs"""
import sys, time; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
shapes = [[int(x) for x in s.split(',')] for s in sys.argv[1].split(';')]
NS, STEPS, MAXS = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
t0 = time.time(); POOL = {}
MINP = {1: 1, 2: 2, 3: 4, 4: 6, 5: 10, 6: 16, 7: 28, 8: 48, 9: 112}
def getpool(d):
    if d not in POOL:
        if d == 9: POOL[d] = ([()], np.ones((1, NP), bool))
        else: POOL[d] = pool(d, minports=MINP[d])
    return POOL[d]
rng = np.random.default_rng(int(time.time() * 1000) % 2**31)
for dims in shapes:
    for d in set(dims): getpool(d)
    bestf = None; budget = MAXS / len(shapes); ts = time.time()
    for seed in range(NS):
        if time.time() - ts > budget: break
        cur = [rng.integers(len(POOL[d][0])) for d in dims]
        S = np.array([POOL[d][1][i] for d, i in zip(dims, cur)])
        f = int(hall_deficit(S).clip(0).sum()); T = 3.0; fb = f
        for step in range(STEPS):
            j = rng.integers(len(dims)); d = dims[j]
            if len(POOL[d][0]) == 1: continue
            new = rng.integers(len(POOL[d][0]))
            S2 = S.copy(); S2[j] = POOL[d][1][new]
            f2 = int(hall_deficit(S2).clip(0).sum())
            if f2 <= f or rng.random() < np.exp((f - f2) / T): cur[j] = new; S = S2; f = f2
            fb = min(fb, f); T = max(0.05, T * 0.9995)
            if f == 0: break
        bestf = fb if bestf is None else min(bestf, fb)
        if f == 0:
            print('  HALL PASSED: dims %s duals %s' % (dims, [POOL[d][0][i] for d, i in zip(dims, cur)]), flush=True); break
    print('shape %s cost %d: least Hall deficit %d (%d seeds x <= %d steps)  %.0fs' % (dims, sum(dims), bestf, seed + 1, STEPS, time.time() - t0), flush=True)
