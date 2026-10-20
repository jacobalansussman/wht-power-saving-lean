"""subunit_onerow.py (d-R7): on U' (35 ports, root rank 7): which single frames F (subspaces of V) can carry one row
of a rank-7 (K=7) fitting decoder whose other 6 rows are unrestricted?  Masked ALS, 2 restarts per F.
usage: subunit_onerow.py dim max_frames K"""
import sys, time, itertools; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
from als import als
d, MAXF, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
U = np.array(lines_of_dual([1, 2])); _O = ORTH[np.ix_(U, U)]; U = U[_O.sum(1) < len(U) - 1]
OU = ORTH[np.ix_(U, U)]; nU = len(U)
GU = GP[np.ix_(U, U)].astype(float); rU = qrank(GU)
Vb = f2basis([x for x in range(512) if par(x & 1) == 0 and par(x & 2) == 0])
vecs = [sum(Vb[i] for i in range(7) if c >> i & 1) for c in range(128)]
frames = []
for piv in itertools.combinations(range(7), d):
    frees = [[j for j in range(p + 1, 7) if j not in piv] for p in piv]
    for bits in itertools.product(*[range(1 << len(f)) for f in frees]):
        rows = [(1 << p) | sum(1 << j for i, j in enumerate(f) if b >> i & 1) for p, f, b in zip(piv, frees, bits)]
        span = {0}
        for r in rows: span |= {x ^ r for x in span}
        labs = {vecs[c] for c in span}
        frames.append(np.array([LAB[t] in labs for t in U]))
rng = np.random.default_rng(d)
idx = rng.permutation(len(frames))[:MAXF]
t0 = time.time(); ok = []; nat_ok = 0; lin = 0
for i in idx:
    F = frames[i]
    # natural possible? roots of U' outside F of rank < rU
    if (~F).any() and qrank(LAM[U][~F]) < rU: lin += 1
    M = np.ones((K, nU), bool); M[-1] = F
    best = min(als(M, OU, iters=600, rng=rng, ridge=1e-13)[0] for _ in range(2))
    # the row must be used: check that removing it breaks the fit is implicit only if K = rank; report residual
    if best < 1e-7: ok.append((int(i), int(F.sum())))
print('dim %d frames tested %d (ports per frame %s): one restricted row fits in %d; frames admitting a linear total: %d  %.0fs'
      % (d, len(idx), sorted(set(int(frames[i].sum()) for i in idx)), len(ok), lin, time.time() - t0))
print('fitting frames (index, ports):', ok[:20])
