"""split_control.py (d-R7): positive control for masked ALS with dim-7 frames.  Coordinate total o_8 (56 ports, an
8-dim frame) is split into two 7-dim frames F_a, F_b covering its 56 ports (natsplit-style, natural, cost 70);
masked ALS on (7 coordinate masks + Lines(F_a) + Lines(F_b)) must reach ~0.  Then the same with the two 7-dim
masks replaced by random Hall-repairing 7-dim frames (cost 70, non-split)."""
import sys, time; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
from als import als
t0 = time.time(); rng = np.random.default_rng(5)
COORD = np.array([[GRAM[o, P[T]] != 0 for T in range(NP)] for o in list(O0)])
C8 = COORD[7]
ds7, M7 = pool(7)
inC = M7 & C8[None, :]
cnt = inC.sum(1)
best = None
cand = np.nonzero(cnt >= 20)[0]
for a in cand:
    rest = C8 & ~M7[a]
    ok = np.nonzero((M7[:, rest].all(1)))[0]
    if len(ok): best = (a, ok[0]); break
a, b = best
print('split of coordinate total 8 into 7-dim frames: ports of C8 in F_a %d, in F_b %d, union covers C8: %s' % (inC[a].sum(), inC[b].sum(), bool((C8 <= (M7[a] | M7[b])).all())))
for name, masks in [('split (natural, masks = Lines(F) cap C8)', np.vstack([COORD[:7], (M7[a] & C8)[None], (M7[b] & C8)[None]])),
                    ('split (masks = full Lines(F_a), Lines(F_b))', np.vstack([COORD[:7], M7[a][None], M7[b][None]]))]:
    for rs in range(3):
        res, X, R, h = als(masks, ORTH, iters=1500, rng=rng)
        print('%s restart %d: residual %.2e nat %s  %.0fs' % (name, rs, res, np.round(natural_test(X @ R), 4), time.time() - t0), flush=True)
