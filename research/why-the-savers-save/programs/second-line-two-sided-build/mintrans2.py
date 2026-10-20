#!/usr/bin/env python3
"""an-build: two-sided slack on a gcert/1 local word WITH the dirt step (own code; exact LP, checked integral).
Frame bit b (s in the frame) per gate; plus a time-0 node v_r for every dirt-step slot (K = Mt[y, r] != 0, from
dirtjoin.py) and a cut node c_k for every retained slot (its frame f when the copy k takes it).
  x gates: b = 0;  y gates: b = 1;  chains monotone along every register: [v_r] -> A gates -> [c_k] -> B gates.
  v_r = 1: the slot meets the y's at time 0 in their frame (holds s from then on);
  v_r = 0: the dirt step of r goes through an erased crossing copy  <s> -> 0  (rank 1): +1 in every GB call
           (in GF calls the slot must instead be clean, i.e. banked, or pay the same +1).
  c_k = 1: the copy's transition f+s <-> <s> is nested (rank dim f, as one-sided); c_k = 0: f <-> <s> non-nested, +1 per call.
Objective (moves per replica): 4 * #(v_r = 0) + 9 * #(c_k = 0), then (tie-break) the number of switching slots.
usage: python3 -I mintrans2.py FILE DIRTJSON [out.json]"""
import gzip, json, sys, os
from collections import Counter
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
d = json.load(gzip.open(sys.argv[1])); dj = json.load(open(sys.argv[2]))
v, R = d['v'], d['R']; n = 2 * v + R
dim = [len(f) for f in d['frames']]
A, B = d['A'], d['B']; GA = A + B; nA = len(A); NG = len(GA)
tot = [s_ for _, s_, _ in d['ret']]
def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
RG = [regs(g) for g in GA]
K = dj['Kslots']
vnode = {r: NG + i for i, r in enumerate(K)}
cnode = {t: NG + len(K) + i for i, t in enumerate(tot)}
NV = NG + len(K) + len(tot)
chain = {r: ([vnode[r]] if r in vnode else []) for r in range(n)}
for i, rs in enumerate(RG):
    if i == nA:
        for t in tot: chain[t].append(cnode[t])
    for r in rs: chain[r].append(i)
lo = np.zeros(NV); hi = np.ones(NV)
for i, rs in enumerate(RG):
    if any(r < v for r in rs): hi[i] = 0
    if any(v <= r < 2 * v for r in rs): lo[i] = 1
assert (lo <= hi).all()
rows, cols, vals = [], [], []; k = 0
for r in range(n):
    c = chain[r]
    for a, b in zip(c, c[1:]):
        rows += [k, k]; cols += [a, b]; vals += [1, -1]; k += 1
Aub = coo_matrix((vals, (rows, cols)), shape=(k, NV)).tocsr()
slots = list(range(2 * v, n))
def solve(tie):
    cost = np.zeros(NV)
    for r in K: cost[vnode[r]] -= 4.0            # minimise 4*(1 - v_r)  <=>  maximise v_r
    for t in tot: cost[cnode[t]] -= 9.0
    if tie:
        for r in slots:
            c = chain[r]
            if len(c) > 1 and c[0] != c[-1]:
                cost[c[-1]] += 1e-3; cost[c[0]] -= 1e-3
    res = linprog(cost, A_ub=Aub, b_ub=np.zeros(k), bounds=list(zip(lo, hi)), method='highs')
    assert res.status == 0, res.message
    x = res.x
    integral = bool(np.all(np.minimum(np.abs(x), np.abs(x - 1)) < 1e-9))
    return np.round(x).astype(int), integral
b, integral = solve(True)
for r in range(n):
    c = chain[r]; assert all(b[p] <= b[q] for p, q in zip(c, c[1:]))
cross = [r for r in K if b[vnode[r]] == 0]
nonnested = [t for t in tot if b[cnode[t]] == 0]
cls = {}
for r in slots:
    c = chain[r]
    if not c: cls[r] = 'free'
    elif b[c[0]] == 1: cls[r] = 'Y'
    elif b[c[-1]] == 0: cls[r] = 'X'
    else: cls[r] = 'S'
# exact replay of the climbs with the s-bit (a switching climb pays +1); copies: dim f, +1 if non-nested
start, final = d['start'], d['final']
cur = list(start); sb = [False] * n
for r in range(v, 2 * v): sb[r] = True
for r in slots:
    if cls[r] == 'Y': sb[r] = True
H = Counter(); Hc = Counter()
def climb(r, f, s_need):
    o = cur[r]; gain = dim[f] - dim[o] + (1 if (s_need and not sb[r]) else 0)
    assert dim[f] >= dim[o] and (sb[r] <= s_need), (r, o, f, sb[r], s_need)
    if gain > 0: H[gain] += 1
    cur[r] = f; sb[r] = sb[r] or s_need
for i, rs in enumerate(RG):
    if i == nA:
        for kk, s_, f in d['ret']:
            climb(s_, f, bool(b[cnode[s_]])); Hc[dim[f] + (0 if b[cnode[s_]] else 1)] += 1
    for r in rs: climb(r, GA[i][1], bool(b[i]))
for r in range(n): climb(r, final[r], sb[r])
span = {r: dim[final[r]] - dim[start[r]] + (1 if cls[r] == 'S' else 0) for r in slots}
C = Counter(cls.values())
wid = {c: Counter(span[r] for r in slots if cls[r] == c) for c in C}
prof = Counter(H); prof.update(Hc)
out = dict(file=sys.argv[1].split('/')[-1], integral_lp=integral, Kslots=len(K), crossing_copies_per_GB_call=len(cross),
           crossing_slots_meeting_x_in_A=len(set(cross) & set(dj['K_meet_x_in_A'])),
           nonnested_copies_per_call=len(nonnested), extra_per_replica=4 * len(cross) + 9 * len(nonnested),
           switching=C.get('S', 0), classes=dict(C), widths={c: dict(sorted(w.items())) for c, w in wid.items()},
           width_sum={c: sum(a * m for a, m in w.items()) for c, w in wid.items()},
           cross_slot_classes=dict(Counter(cls[r] for r in cross)),
           profile_GF=dict(sorted(prof.items())), rank_GF=sum(a * m for a, m in prof.items()), cst=sum(a * m for a, m in Hc.items()))
out['rank_minus_widths'] = out['rank_GF'] - sum(out['width_sum'].values())
out['D_per_replica_upper'] = 2 * 896 - 9 * 64 - out['extra_per_replica']
print(json.dumps(out))
if len(sys.argv) > 3: json.dump(out, open(sys.argv[3], 'w'), indent=1)
