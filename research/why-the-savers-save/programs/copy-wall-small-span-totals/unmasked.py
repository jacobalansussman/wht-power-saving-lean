"""unmasked.py (d-R7): rank-K fitting matrices on P112 with no read masks, random restarts; is each one natural
(rowspace contains Lambda D)?  usage: unmasked.py K restarts iters"""
import sys, time; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import warnings; warnings.filterwarnings('ignore')
from model import *
from als import als
import numpy as np
K, NR, IT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(K * 1000 + 7)
t0 = time.time()
for rs in range(NR):
    res, X, R, h = als(np.ones((K, NP), dtype=bool), ORTH, iters=IT, rng=rng, ridge=1e-14)
    B = X @ R
    sv = np.linalg.svd(B, compute_uv=False)
    print('K=%d rs %d res %.2e  sv[7..] %s  natural(r, s_min, s_next, min|d|) %s  %.0fs' % (K, rs, res, np.round(sv[7:K+1] / sv[0], 5), np.round(natural_test(B), 6), time.time() - t0), flush=True)
