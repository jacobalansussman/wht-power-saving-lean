"""famsearch.py (d-R7): local search for K frames (dims given) whose read masks pass octad Hall on P112.
usage: famsearch.py dims(comma) seeds steps out.npz"""
import sys, time; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
dims = [int(x) for x in sys.argv[1].split(',')]; NS = int(sys.argv[2]); STEPS = int(sys.argv[3]); outp = sys.argv[4]
t0 = time.time()
POOL = {d: pool(d, minports=int(__import__("os").environ.get("MINP_%d" % d, {6: 16, 7: 28, 8: 48}[d]))) for d in set(dims)}
print('pools', {d: len(POOL[d][0]) for d in POOL}, '%.0fs' % (time.time() - t0), flush=True)
rng = np.random.default_rng(int(time.time()) % 100000)
found = []
for seed in range(NS):
    cur = [rng.integers(len(POOL[d][0])) for d in dims]
    S = np.array([POOL[d][1][i] for d, i in zip(dims, cur)])
    f = int(hall_deficit(S).clip(0).sum()); best = f
    T = 3.0
    for step in range(STEPS):
        j = rng.integers(len(dims)); d = dims[j]
        new = rng.integers(len(POOL[d][0]))
        S2 = S.copy(); S2[j] = POOL[d][1][new]
        f2 = int(hall_deficit(S2).clip(0).sum())
        if f2 <= f or rng.random() < np.exp((f - f2) / T):
            cur[j] = new; S = S2; f = f2
        T = max(0.05, T * 0.9995)
        if f == 0: break
    print('seed %d dims %s: Hall failures %d after %d steps, ports read per total %s, uncovered %d  %.0fs' % (
        seed, dims, f, step + 1, S.sum(1).tolist(), int((S.sum(0) == 0).sum()), time.time() - t0), flush=True)
    if f == 0: found.append([POOL[d][0][i] for d, i in zip(dims, cur)])
np.save(outp, np.array(found, dtype=object), allow_pickle=True)
print('families passing Hall:', len(found))
