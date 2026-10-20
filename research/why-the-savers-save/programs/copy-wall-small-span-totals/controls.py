"""controls.py (d-R7): ALS controls on P112.  (a) the 8 coordinate-total masks (must fit, residual -> 0);
(b) unmasked rank 8 and rank 9 (natural test on what is found)."""
import sys, time; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import warnings; warnings.filterwarnings('ignore')
from model import *
from als import als
import numpy as np
t0 = time.time()
# coordinate totals: the removed octad lines o; total reads ports non-orthogonal to o
Oidx = list(O0)
RM8 = np.array([[GRAM[o, P[T]] != 0 for T in range(NP)] for o in Oidx])
print('coordinate masks: reads per total', RM8.sum(1), 'frame dims', [f2dim([LAB[T] for T in range(NP) if RM8[k, T]]) for k in range(8)])
rng = np.random.default_rng(1)
for rs in range(3):
    res, X, R, h = als(RM8, ORTH, iters=400, rng=rng)
    print('(a) coord masks restart', rs, 'residual %.2e' % res, 'natural:', natural_test(X @ R)[:3], round(time.time() - t0, 1), flush=True)
full = lambda K: np.ones((K, NP), dtype=bool)
for K in (8, 9, 10):
    for rs in range(int(sys.argv[1]) if len(sys.argv) > 1 else 3):
        res, X, R, h = als(full(K), ORTH, iters=600, rng=rng)
        B = X @ R
        print('(b) unmasked K=%d restart %d residual %.2e rank/natdim/min|d| %s  hist %s' % (K, rs, res, natural_test(B)[:3], ['%.1e' % x for x in h[::6]]), round(time.time() - t0, 1), flush=True)
