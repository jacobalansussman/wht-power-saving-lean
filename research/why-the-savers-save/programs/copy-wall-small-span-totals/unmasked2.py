"""unmasked2.py (d-R7): many restarts of unmasked rank-K ALS on P112 (several init styles); classify each fitting
B: natural (Lambda D in rowspace with all d_T != 0), partially natural (some d_T ~ 0), or non-natural.
usage: unmasked2.py K restarts iters [init: gauss|sparse|mixed]"""
import sys, time; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import warnings; warnings.filterwarnings('ignore')
from model import *
from als import als
import numpy as np, collections
K, NR, IT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); init = sys.argv[4] if len(sys.argv) > 4 else 'gauss'
rng = np.random.default_rng(K * 1000 + 11 + hash(init) % 1000)
t0 = time.time(); cls = collections.Counter(); worst = []
for rs in range(NR):
    R0 = rng.standard_normal((K, NP))
    if init == 'sparse': R0 *= rng.random((K, NP)) < 0.3
    if init == 'mixed': R0[:8] = (LAM @ rng.standard_normal((8, 8))).T; R0[8:] *= 5
    res, X, R, h = als(np.ones((K, NP), dtype=bool), ORTH, iters=IT, rng=rng, ridge=1e-14, R0=R0)
    if res > 1e-8: cls['not converged'] += 1; worst.append(res); continue
    r, smin, snext, dmin = natural_test(X @ R)
    c = 'natural' if smin < 1e-6 and dmin > 1e-3 else ('partial' if smin < 1e-6 else 'NON-NATURAL')
    cls[(c, r)] += 1
    if c != 'natural': print('rs %d res %.1e rank %d smin %.2e snext %.2e dmin %.2e' % (rs, res, r, smin, snext, dmin), flush=True)
print('K=%d init=%s restarts=%d: %s  unconverged residuals %s  %.0fs' % (K, init, NR, dict(cls), ['%.1e' % x for x in sorted(worst)[:5]], time.time() - t0))
