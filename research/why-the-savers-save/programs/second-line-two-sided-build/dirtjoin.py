#!/usr/bin/env python3
"""an-build: the dirt step of the Lean word (Work/Carrier/Live.lean g1 M = 1 - PY M PS; Phased.lean f1 = MA . g1(Mt)):
"y -= K s" at the very start, where K = Mt restricted to (y rows, slot columns).  Its slots meet the y's at frame 0 at
time 0 (the erased invocation starts every slot at zeroL).  In the backward word its transpose  s += K^T y  is the ONLY
entry of y information before phase A.  Exact replay (Fractions) of the forward map Mt of a gcert/1 word: A, the scatter
(retained slot k copied into total k, y += Jr * totals), B.  Reports which slots have K != 0, whether they meet an x in
phase A, and checks the Lean facts hx, hy, hid (x rows unit, y columns unit, y_S <- x_t = [S = t]).  Also checks that the
backward word built from these pieces delivers y information to every x.
usage: python3 -I dirtjoin.py FILE"""
import gzip, json, sys
from fractions import Fraction as F
d = json.load(gzip.open(sys.argv[1]))
v, R = d['v'], d['R']; n = 2 * v + R
A, B = d['A'], d['B']
def adds(g):
    return [(t, g[2], F(a, b)) for t, a, b in g[3]] if g[0] == 'out' else [(g[2], s, F(a, b)) for s, a, b in g[3]]
row = [{r: F(1)} for r in range(n)]          # content of register r as a combination of inputs
def addto(t, s, c):
    for k, val in row[s].items():
        nv = row[t].get(k, 0) + c * val
        if nv: row[t][k] = nv
        else: row[t].pop(k, None)
for g in A:
    for t, s, c in adds(g): addto(t, s, c)
tot = {k: s for k, s, _ in d['ret']}
copies = {k: dict(row[s]) for k, s in tot.items()}
for T, entries in enumerate(d['scat']['table']):
    for k, a, b in entries:
        c = F(a, b)
        for key, val in copies[k].items():
            nv = row[v + T].get(key, 0) + c * val
            if nv: row[v + T][key] = nv
            else: row[v + T].pop(key, None)
for g in B:
    for t, s, c in adds(g): addto(t, s, c)
isx = lambda r: r < v; isy = lambda r: v <= r < 2 * v; iss = lambda r: r >= 2 * v
hx = all(row[t] == {t: 1} for t in range(v))
hy = all(all(not (k == y and r != y) for r in range(n) for k in row[r]) for y in range(v, 2 * v))
hid = all(row[v + S].get(t, 0) == (1 if S == t else 0) for S in range(v) for t in range(v))
print('Lean facts on the replayed forward map: hx', hx, '| y columns unit', hy, '| hid', hid)
Kslots = sorted({k for T in range(v) for k in row[v + T] if iss(k)})
print('slots with K = Mt[y, slot] != 0 (joined with the y''s at time 0 by the dirt step):', len(Kslots), 'of', R)
xA = set()
for g in A:
    rs = [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
    if any(isx(r) for r in rs): xA.update(r for r in rs if iss(r))
print('  of them meeting an x in phase A directly:', len(set(Kslots) & xA))
# backward word = h8^T . MB^-T . scat^-T . MA^-T . g1i^T ; follow y information
has = [False] * n
for r in range(v, 2 * v): has[r] = True
for k in Kslots: has[k] = True                     # g1i^T: s += K^T y
for g in A:
    for t, s, c in adds(g):
        if has[t]: has[s] = True
yinfo_cut = any(has[r] for r in range(v, 2 * v))
for k, s in tot.items(): has[s] = has[s] or yinfo_cut
for g in B:
    for t, s, c in adds(g):
        if has[t]: has[s] = True
print('backward word with the dirt step: x reached by y information:', sum(has[r] for r in range(v)), 'of', v)
json.dump(dict(Kslots=Kslots, K_meet_x_in_A=sorted(set(Kslots) & xA), hx=hx, hid=hid), open('dirtjoin_' + sys.argv[1].split('/')[-1].split('.')[0] + '.json', 'w'))
