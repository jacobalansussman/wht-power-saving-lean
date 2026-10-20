# fulltouch3.py -- fulltouch.py's sequences on ONE two-row array in a NON-ABELIAN letter group (order 8, two involutions)
import numpy as np, sys
from fulltouch import is_inv, close, idx, run
from fam import POOL
invs = [u for u in POOL if is_inv(u)]
rng = np.random.default_rng(int(sys.argv[1]))
found = 0
for trial in range(2000):
    u1, u2 = invs[int(rng.integers(len(invs)))], invs[int(rng.integers(len(invs)))]
    if np.allclose(u1 @ u2, u2 @ u1): continue
    els = close([u1, u2], cap=8)
    if els is None or len(els) != 8: continue
    i0, i1, i2 = 0, idx(els, u1), idx(els, u2); i12 = idx(els, u1 @ u2)
    seqs = {
      '(1+u1), then 1': [{i0: 1, i1: 1}, {i0: 1}],
      '(1+u1), (1+u2), then 1': [{i0: 1, i1: 1}, {i0: 1, i2: 1}, {i0: 1}],
      '(1+u1)(1+u2), then (1-u1), then 1': [{i0: 1, i1: 1, i2: 1, i12: 1}, {i0: 1, i1: -1}, {i0: 1}],
      'projector of Gamma, then 1': [{i: 1.0 / len(els) for i in range(len(els))}, {i0: 1}],
      '(1+u1+u2), then 1': [{i0: 1, i1: 1, i2: 1}, {i0: 1}],
    }
    for lab, sq in seqs.items(): run(els, sq, 'D8 ' + lab)
    found += 1
    if found >= 1: break
print('DONE')
