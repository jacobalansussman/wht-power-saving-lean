#!/usr/bin/env python3
"""vpairs_sanity.py: (a) the float32 matmul used in vpairs.py equals the exact integer count on a sample of multisets;
(b) positive control for the deletion test: plain 2 in 8 (directions e1..e4 and the columns of h, T = h h^T) passes for
every target; (c) the 7-move enumeration logic refuses 1 in 3 for I (negative control) and accepts 1 in 4 plain."""
import sys, itertools, numpy as np
sys.path.insert(0, '.')
from vcore import setup, pairing, key, targets_and_orth
D = setup(4); P = D['P']; Qs = D['Qs']; GL = D['GL']
T, O = targets_and_orth(D); I = np.eye(4, dtype=np.int64); nI = pairing(Qs, I)
rng = np.random.default_rng(5)
ms = np.array([np.bincount(rng.integers(0, 15, 7), minlength=15) for _ in range(3000)])
a = np.rint(ms.astype(np.float32) @ P.T.astype(np.float32)).astype(np.int64); b = ms @ P.T
print('(a) float32 count == exact count on 3000 random 7-multisets:', np.array_equal(a, b))
def idx(v): return int(sum(int(v[i]) << i for i in range(4))) - 1
ok = 0
for h in GL[rng.choice(len(GL), 300, replace=False)]:
    S = (h @ h.T) % 2
    c = np.zeros(15, dtype=np.int64)
    for i in range(4): c[idx(np.eye(4, dtype=np.int64)[:, i])] += 1; c[idx(h[:, i])] += 1
    ok += bool(np.all(P @ c >= nI + pairing(Qs, S)))
print('(b) plain 2 in 8 passes deletion for %d of 300 random h' % ok)
c = np.zeros(15, dtype=np.int64); c[[0, 1, 3]] = 1
print('(c) 1 in 3 (e1,e2,e3) for I passes deletion?', bool(np.all(P @ c >= nI)), '; 1 in 4 plain passes?', end=' ')
c[7] = 1; print(bool(np.all(P @ c >= nI)))
