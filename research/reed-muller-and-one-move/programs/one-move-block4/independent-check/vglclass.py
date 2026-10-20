#!/usr/bin/env python3
"""vglclass.py (v-t2-block4): merge the O_4-classes of {I, T2, T3} into classes of target multisets under the full
relabelling {T_l} -> {g T_l g^T} (g in GL(4,2)), which leaves the multi-target game (and deletion) unchanged.
Consistency check of the deletion ILP: classes merged together must have equal deletion minima (vtri_low.json has
every class with minimum <= 10)."""
import sys, itertools, json, numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from vcore import setup, key, targets_and_orth
D = setup(4); GL = D['GL']; T, O = targets_and_orth(D); I = np.eye(4, dtype=np.int64)
keys = sorted(T); kidx = {k: i for i, k in enumerate(keys)}
Iidx = kidx[key(I)]
act = np.array([[kidx[key((g @ T[k] @ g.T) % 2)] for k in keys] for g in O])
def o4canon(a, b):
    return min(tuple(sorted((int(act[g, a]), int(act[g, b])))) for g in range(len(O)))
# for each target B pick h with h h^T = B, and g = h^-1
inv = {}
for h in GL:
    k = key((h @ h.T) % 2)
    if k not in inv:
        hi = next(x for x in GL if np.array_equal((x @ h) % 2, I)); inv[k] = hi
def glcanon(trip):   # trip: 3 target indices
    best = None
    for i in range(3):
        g = inv[keys[trip[i]]]
        rest = [kidx[key((g @ T[keys[trip[j]]] @ g.T) % 2)] for j in range(3) if j != i]
        c = o4canon(*rest)
        best = c if best is None or c < best else best
    return best
seen = set(); o4 = []
for a, b in itertools.combinations_with_replacement(range(len(keys)), 2):
    if (a, b) in seen: continue
    orb = {tuple(sorted((int(act[g, a]), int(act[g, b])))) for g in range(len(O))}
    seen |= orb; o4.append(min(orb))
print('O_4-classes of {I,T2,T3}:', len(o4))
grp = defaultdict(list)
for c in o4: grp[glcanon((Iidx, c[0], c[1]))].append(c)
print('classes of target triples under the full relabelling by GL(4,2):', len(grp))
low = json.load(open('vtri_low.json'))
def idx_of(v): return kidx[tuple(v)]
lowM = {o4canon(idx_of(r['T2']), idx_of(r['T3'])): r['M'] for r in low}
bad = 0; dist = Counter()
for G, mem in grp.items():
    ms = {lowM.get(c, '>10') for c in mem}
    if len(ms) > 1: bad += 1
    dist[ms.pop() if len(ms) == 1 else 'MIXED'] += 1
print('GL-classes with mixed deletion minima (should be 0):', bad)
print('GL-class distribution of deletion minimum (<=10 resolved):', sorted(dist.items(), key=str))
# deletion minima of one representative per GL-class (ILP), to split the '>10' group
from scipy.optimize import milp, LinearConstraint, Bounds
from vcore import pairing
P = D['P']; Qs = D['Qs']; nI = pairing(Qs, I)
d2 = Counter()
for G, mem in grp.items():
    a, b = mem[0]
    n = nI + pairing(Qs, T[keys[a]]) + pairing(Qs, T[keys[b]])
    r = milp(c=np.ones(15), constraints=[LinearConstraint(P, lb=n, ub=np.inf)], integrality=np.ones(15),
             bounds=Bounds(0, np.inf), options=dict(mip_rel_gap=0))
    c = np.round(r.x).astype(np.int64); assert r.status == 0 and np.all(P @ c >= n); d2[int(c.sum())] += 1
print('GL-class deletion minima (ILP on one representative each):', sorted(d2.items()))
