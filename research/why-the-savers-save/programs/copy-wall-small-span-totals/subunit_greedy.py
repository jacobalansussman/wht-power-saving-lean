"""subunit_greedy.py (d-R7): on U' (35 ports, root rank 7, natural floor 42 = 7 x 6): start from a 5-dim frame F1 that
carries one row of a rank-7 fitting decoder (subunit_onerow.py), then add 6-dim frames greedily while the masked
ALS still fits (others free).  Goal shape (5, 6^6) = 41 < 42.  usage: subunit_greedy.py seed"""
import sys, time, itertools; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
from als import als
SEED = int(sys.argv[1]); K = 7
U = np.array(lines_of_dual([1, 2])); _O = ORTH[np.ix_(U, U)]; U = U[_O.sum(1) < len(U) - 1]
OU = ORTH[np.ix_(U, U)]; nU = len(U)
Vb = f2basis([x for x in range(512) if par(x & 1) == 0 and par(x & 2) == 0])
vecs = [sum(Vb[i] for i in range(7) if c >> i & 1) for c in range(128)]
def frames_of(d):
    out = []
    for piv in itertools.combinations(range(7), d):
        frees = [[j for j in range(p + 1, 7) if j not in piv] for p in piv]
        for bits in itertools.product(*[range(1 << len(f)) for f in frees]):
            rows = [(1 << p) | sum(1 << j for i, j in enumerate(f) if b >> i & 1) for p, f, b in zip(piv, frees, bits)]
            span = {0}
            for r in rows: span |= {x ^ r for x in span}
            labs = {vecs[c] for c in span}
            out.append(np.array([LAB[t] in labs for t in U]))
    return out
F5 = frames_of(5); F6 = frames_of(6)
rng = np.random.default_rng(SEED); t0 = time.time()
def fits(masks, tries=3):
    M = np.ones((K, nU), bool); M[:len(masks)] = np.array(masks)
    best = (1e9,)
    for _ in range(tries):
        r = als(M, OU, iters=700, rng=rng, ridge=1e-13)
        if r[0] < best[0]: best = r
        if best[0] < 1e-8: break
    return best
start = [1728, 2359, 2224, 2054, 1979, 2380, 2543, 373, 2264, 2179]
for s in start[SEED % len(start):] + start[:SEED % len(start)]:
    chosen = [F5[s]]; dims = [5]
    r = fits(chosen)
    print('start F5[%d] (%d ports): residual %.1e' % (s, F5[s].sum(), r[0]), flush=True)
    if r[0] > 1e-8: continue
    while len(chosen) < K:
        order = rng.permutation(len(F6)); added = False
        for j in order:
            if F6[j].sum() == 0: continue
            r = fits(chosen + [F6[j]], tries=2)
            if r[0] < 1e-8:
                chosen.append(F6[j]); dims.append(6); added = True
                print('   + F6[%d] (%d ports): fits, residual %.1e, dims %s cost so far %d  %.0fs' % (j, F6[j].sum(), r[0], dims, sum(dims), time.time() - t0), flush=True)
                break
        if not added:
            print('   no 6-dim frame extends %s  %.0fs' % (dims, time.time() - t0), flush=True); break
    if len(chosen) == K:
        res, X, R, h = fits(chosen, tries=4)
        print('FULL masked decoder dims %s cost %d residual %.2e' % (dims, sum(dims), res), flush=True)
        np.savez(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/subunit_decoder_%d.npz') % SEED, X=X, R=R, masks=np.array(chosen), U=U)
        break
    if time.time() - t0 > 100: print('time box'); break
