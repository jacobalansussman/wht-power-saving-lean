#!/usr/bin/env python3
"""vtri.py (v-t2-block4): my own classification of target triples {I, T2, T3} (T2, T3 any targets h h^T, repetition
allowed) up to O_4 acting by T -> g T g^T, and the exact integer DELETION minimum of each class (ILP, HiGHS with
mip_rel_gap = 0; every witness re-checked in integers).  Also the LP relaxation of deletion.
Output: distribution of minima; list of classes with minimum <= 10 (with witness counts)."""
import sys, itertools, time, json
import numpy as np
from collections import Counter
from scipy.optimize import milp, LinearConstraint, Bounds
sys.path.insert(0, '.')
from vcore import setup, pairing, key, targets_and_orth
t0 = time.time()
D = setup(4); P = D['P']; Qs = D['Qs']
T, O = targets_and_orth(D)
I = np.eye(4, dtype=np.int64)
keys = sorted(T); kidx = {k: i for i, k in enumerate(keys)}
# action of O_4 on target indices
act = np.array([[kidx[key((g @ T[k] @ g.T) % 2)] for k in keys] for g in O])   # (48, 420)
seen = set(); classes = []
for a, b in itertools.combinations_with_replacement(range(len(keys)), 2):
    if (a, b) in seen: continue
    orb = {tuple(sorted((int(act[g, a]), int(act[g, b])))) for g in range(len(O))}
    seen |= orb; classes.append(((a, b), len(orb)))
print('triple classes {I, T2, T3} up to O_4: %d (orbit sizes sum %d = C(421,2) = %d)' % (
    len(classes), sum(s for _, s in classes), 420 * 421 // 2), flush=True)
nI = pairing(Qs, I); Tn = np.array([pairing(Qs, T[k]) for k in keys])
dist = Counter(); low = []
for (a, b), sz in classes:
    n = nI + Tn[a] + Tn[b]
    res = milp(c=np.ones(15), constraints=[LinearConstraint(P, lb=n, ub=np.inf)], integrality=np.ones(15),
               bounds=Bounds(0, np.inf), options=dict(mip_rel_gap=0))
    assert res.status == 0
    c = np.round(res.x).astype(np.int64); assert np.all(P @ c >= n)
    M = int(c.sum()); dist[M] += 1
    if M <= 10:
        low.append(dict(M=M, T2=[int(v) for v in T[keys[a]].ravel()], T3=[int(v) for v in T[keys[b]].ravel()],
                        size=sz, counts=c.tolist()))
print('deletion minimum distribution:', sorted(dist.items()))
print('classes allowing 11 moves (minimum <= 11):', sum(v for k, v in dist.items() if k <= 11))
json.dump(low, open('vtri_low.json', 'w'))
print('time %.1f s' % (time.time() - t0))
