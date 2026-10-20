#!/usr/bin/env python3
"""an-build: the fewest s-switching slots of a gcert/1 local word under two-sided slack, exactly (own code).
One bit b_i per gate (s in the common frame of gate i).  Fixed: b = 0 on every gate naming an x (x ends in s^perp),
b = 1 on every gate naming a y (y starts on <s>); every total holds s at the cut (virtual node fixed to 1 between its
last A gate and its first B gate).  Monotone along every register: b_i <= b_j for consecutive gates i < j of it.
A slot SWITCHES iff b(first gate) = 0 and b(last gate) = 1; minimise  sum_slots b_last - b_first.
The constraint matrix is a difference system (totally unimodular), so the LP optimum is integral (checked).
Then: exact replay of the frames with the s-bit -> local block profile (a switching climb pays +1), residual widths by
host class (X: never s; Y: s from the start; S: switches once).  usage: python3 -I mintrans.py FILE [out.json]"""
import gzip, json, sys, os
from collections import Counter
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict(a.lstrip('-').split('=') for a in sys.argv[1:] if a.startswith('--'))
CUT = opts.get('cut', 'totals')      # 'totals': ret slots hold s from the cut (first model); 'copies': ret slots s-free at the cut, the C copies carry s (non-nested step, +1 each)
TIE = opts.get('tie', 'none')        # 'x': among optimal solutions prefer s-free gates (more X-class width); 'y': prefer s
d = json.load(gzip.open(args[0]))
h, v, R = d['h'], d['v'], d['R']; n = 2 * v + R
dim = [len(f) for f in d['frames']]
A, B = d['A'], d['B']; GA = A + B; nA = len(A); NG = len(GA)
tot = [s_ for _, s_, _ in d['ret']]
def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
RG = [regs(g) for g in GA]
# variables: gates 0..NG-1, then one cut node per total
cutvar = {t: NG + k for k, t in enumerate(tot)}
NV = NG + len(tot)
chain = {r: [] for r in range(n)}
for i, rs in enumerate(RG):
    if i == nA:
        for t in tot: chain[t].append(cutvar[t])
    for r in rs: chain[r].append(i)
if nA == NG:
    for t in tot: chain[t].append(cutvar[t])
lo = np.zeros(NV); hi = np.ones(NV)
for t in tot:
    if CUT == 'totals': lo[cutvar[t]] = 1
    else: hi[cutvar[t]] = 0
for i, rs in enumerate(RG):
    if any(r < v for r in rs): hi[i] = 0
    if any(v <= r < 2 * v for r in rs): lo[i] = 1
assert (lo <= hi).all(), 'a gate names both an x and a y (or an x after a total joined it): two-sided is illegal at this order'
rows, cols, vals = [], [], []; k = 0
for r in range(n):
    c = chain[r]
    for a, b in zip(c, c[1:]):
        if a == b: continue
        rows += [k, k]; cols += [a, b]; vals += [1, -1]; k += 1            # b_a - b_b <= 0
Aub = coo_matrix((vals, (rows, cols)), shape=(k, NV)).tocsr(); bub = np.zeros(k)
cost = np.zeros(NV)
slots = list(range(2 * v, n))
for r in slots:
    c = chain[r]
    if c and c[0] != c[-1]:
        cost[c[-1]] += 1; cost[c[0]] -= 1
BIG = 10.0 ** 6
if TIE != 'none':
    wgt = np.zeros(NV)
    for r in slots:
        c = chain[r]
        if c: wgt[c[-1]] += (dim[d['final'][r]] - dim[d['start'][r]])
    cost = BIG * cost + (wgt if TIE == 'x' else -wgt)
res = linprog(cost, A_ub=Aub, b_ub=bub, bounds=list(zip(lo, hi)), method='highs')
assert res.status == 0, res.message
from collections import Counter as _C
sw_check = None
b = res.x
integral = bool(np.all(np.minimum(np.abs(b), np.abs(b - 1)) < 1e-9))
b = np.round(b).astype(int)
# feasibility re-check of the rounded solution
for r in range(n):
    c = chain[r]
    assert all(b[x] <= b[y] for x, y in zip(c, c[1:]))
cls = {}
for r in slots:
    c = chain[r]
    if not c: cls[r] = 'free'
    elif b[c[0]] == 1: cls[r] = 'Y'
    elif b[c[-1]] == 0: cls[r] = 'X'
    else: cls[r] = 'S'
# exact replay with the s-bit: frames are (frame id, s); a climb pays dim difference + [s gained]
start, final = d['start'], d['final']
cur = list(start); sb = [False] * n
for r in range(v, 2 * v): sb[r] = True
for r in slots:
    if cls[r] == 'Y': sb[r] = True
H = Counter(); Hc = Counter()
def climb(r, f, s_need, kind='g'):
    o = cur[r]; gain = dim[f] - dim[o] + (1 if (s_need and not sb[r]) else 0)
    assert dim[f] >= dim[o] and (sb[r] <= s_need or kind == 'final'), (r, o, f)
    if gain > 0: H[gain] += 1
    cur[r] = f; sb[r] = sb[r] or s_need
for i, rs in enumerate(RG):
    if i == nA:
        for _, s_, f in d['ret']:
            if CUT == 'totals':
                climb(s_, f, True); Hc[dim[f]] += 1        # total holds f + s; descends to <s>; billed dim f
            else:
                climb(s_, f, False); Hc[dim[f] + 1] += 1   # copy: f (no s) <-> <s>, non-nested, rank dim f + 1
    for r in rs: climb(r, GA[i][1], bool(b[i]))
if nA == NG:
    for _, s_, f in d['ret']: climb(s_, f, True); Hc[dim[f]] += 1
span = {}
for r in range(n):
    s_end = sb[r] if r >= 2 * v else (r >= v)            # slots end at their final frame keeping their s-bit
    o_dim = dim[start[r]] + (1 if (r >= v and r < 2 * v) else 0) + (1 if (r >= 2 * v and cls[r] == 'Y') else 0)
    climb(r, final[r], sb[r], 'final')
for r in slots:
    span[r] = dim[final[r]] - dim[start[r]] + (1 if cls[r] == 'S' else 0)
prof = Counter(H); prof.update(Hc)
C = Counter(cls.values())
wid = {c: Counter(span[r] for r in slots if cls[r] == c) for c in C}
out = dict(file=args[0].split('/')[-1], cut=CUT, tie=TIE, v=v, R=R, integral_lp=integral, lp_objective=round(float(np.round(res.x) @ np.array([(1 if i in set() else 0) for i in range(NV)])),3) if False else None, switching=C.get('S', 0),
           classes=dict(C), widths={c: dict(sorted(w.items())) for c, w in wid.items()},
           width_sum={c: sum(a * m for a, m in w.items()) for c, w in wid.items()},
           totals=dict(Counter(cls[t] for t in tot)), profile=dict(sorted(prof.items())), rank=sum(a * m for a, m in prof.items()),
           cst=sum(a * m for a, m in Hc.items()))
out['rank_minus_widths'] = out['rank'] - sum(out['width_sum'].values())
print(json.dumps(out))
if len(args) > 1: json.dump(out, open(args[1], 'w'), indent=1)
