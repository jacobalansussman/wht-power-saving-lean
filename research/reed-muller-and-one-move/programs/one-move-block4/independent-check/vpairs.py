#!/usr/bin/env python3
"""vpairs.py (v-t2-block4): my own check of the pair claims at block 4.
(1) number of targets h h^T, |O_4|, number of O_4-classes of T (pairs {I, T});
(2) for EVERY one of the 420 targets T: integer deletion minimum of {I, T} by ILP (HiGHS), witness re-checked exactly;
(3) independent of any solver: EXHAUSTIVE enumeration of all multisets of 7 directions (C(21,7) = 116280); for each,
    check deletion for {I, T} for every T. Claim: none passes (so no 2 in 7 by deletion alone, for any T).
Deletion (multi-target slice): for every nonzero form Q, #moves z with Q(z) = 1 >= #arrays l with <Q, T_l> = 1."""
import sys, itertools, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
sys.path.insert(0, '.')
from vcore import setup, pairing, key, targets_and_orth
t0 = time.time()
D = setup(4); P = D['P']; Qs = D['Qs']
T, O = targets_and_orth(D)
I = np.eye(4, dtype=np.int64)
print('|GL| = %d, targets = %d, |O_4| = %d' % (len(D['GL']), len(T), len(O)))
# classes of T under O_4
seen = set(); classes = []
for k, S in T.items():
    if k in seen: continue
    orb = {key((g @ S @ g.T) % 2) for g in O}
    seen |= orb; classes.append((k, len(orb)))
print('O_4-classes of T (pairs {I,T}): %d; sizes %s' % (len(classes), sorted(sz for _, sz in classes)))
nI = pairing(Qs, I)
def delmin(n):
    res = milp(c=np.ones(15), constraints=[LinearConstraint(P, lb=n, ub=np.inf)], integrality=np.ones(15),
               bounds=Bounds(0, np.inf), options=dict(mip_rel_gap=0))
    c = np.round(res.x).astype(np.int64)
    assert res.status == 0 and np.all(P @ c >= n)
    return int(c.sum()), c
dm = {}
for k, S in T.items():
    n = nI + pairing(Qs, S)
    M, c = delmin(n); dm[k] = M
from collections import Counter
print('ILP deletion minimum over all 420 targets T (pair {I,T}):', sorted(Counter(dm.values()).items()))
print('I alone: deletion minimum', delmin(nI)[0], ' {I,I}:', delmin(2 * nI)[0])
# exhaustive multisets of 7 directions
ms = np.array([np.bincount(np.array(cmb), minlength=15) for cmb in itertools.combinations_with_replacement(range(15), 7)])
print('multisets of 7 directions:', len(ms))
Pf = P.T.astype(np.float32)
Tn = np.array([pairing(Qs, S) for S in T.values()]).astype(np.int8)  # (420, 1023)
need = nI[None, :].astype(np.int8) + Tn                        # (420, 1023)
nokI = 0; npass = 0
for a in range(0, len(ms), 5000):
    cnt = np.rint(ms[a:a + 5000].astype(np.float32) @ Pf).astype(np.int8)   # |D_Q|
    okI = np.all(cnt >= nI[None, :], axis=1)
    nokI += int(okI.sum())
    for row in cnt[okI]:
        npass += int(np.all(row[None, :] >= need, axis=1).sum())
print('multisets passing deletion for I alone:', nokI)
print('(multiset, T) pairs passing deletion for {I, T} with 7 moves: %d  (claim: 0)' % npass)
# also 8 moves sanity: count T for which some 8-multiset passes (should be all 420 since ILP gives 8)
print('time %.1f s' % (time.time() - t0))
