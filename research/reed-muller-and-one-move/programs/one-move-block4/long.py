#!/usr/bin/env python3
"""long.py (t2-block4): one long LM run of the multi-target slice (lms.py) on a FIXED order of directions; prints cost
and |x| every 50 iterations (border behaviour: cost -> 0 with |x| -> infinity would be an approximate unit).
Usage: long.py m targets dirs seed maxtime [x0.npy]"""
import sys, time, numpy as np
sys.path.insert(0, '.')
import lms
np.seterr(all='ignore')
m = int(sys.argv[1]); T = [lms.parse_T(t, m) for t in sys.argv[2].split(',')]; dirs = [int(v) for v in sys.argv[3].split(',')]
seed = int(sys.argv[4]); maxtime = float(sys.argv[5])
S = lms.Slice(m, T, dirs); rng = np.random.default_rng(seed)
x = np.load(sys.argv[6]) if len(sys.argv) > 6 else (rng.standard_normal(S.nu) + 1j * rng.standard_normal(S.nu)) * 0.5
t0 = time.time(); k = 0
while time.time() - t0 < maxtime:
    x, cost, it = lms.lm(S, x, iters=50, stall=0.0, maxtime=maxtime)
    k += 50
    print('it %5d cost %.8f |x| %.3e max|x| %.3e t %.0f' % (k, cost, np.linalg.norm(x), np.abs(x).max(), time.time() - t0), flush=True)
    if cost < 1e-20: break
np.save('long_%d.npy' % seed, x)
