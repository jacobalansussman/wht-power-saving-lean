#!/usr/bin/env python3
"""u-helpers-hub: what the 112 U8 helpers that never visit a packet plane do (entering helpers, start dims 1/4/5).
Replays the CANONICAL word over Q from symbolic sources (helpers start at 0 = clean columns) and reports, for each
such helper: its entrance frame (and whether it is a source label <u_s>), the phases that write it / read it,
its value when first read (as a source combination), and its readers' classes.
Usage: python3 -I u8_nonplane.py CANONICAL-WORD.json.gz"""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
B = [basis(f) for f in F]; dims = [len(b) for b in B]
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
start = c['start']
ent = [r for r in range(2*v, n) if dims[start[r]] > 0]
# plane visitors (from the census rule): recompute quickly
cur = list(start); vis3 = Counter()
def inside(a, b):
    for x in B[a]:
        y = x
        for z in B[b]: y = min(y, y ^ z)
        if y: return False
    return True
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not inside(cur[r], f): return
    cur[r] = f
    if dims[f] == 3: vis3[r] += 1
order = ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']
for ph in order:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        climb(e[1], e[-1])
        if e[0] == 'add': climb(e[2], e[-1])
nonplane = [r for r in range(2*v, n) if not vis3[r]]
print('non-plane helpers', len(nonplane), 'entering among them', sum(1 for r in nonplane if r in ent))
print('start dims of non-plane helpers', Counter(dims[start[r]] for r in nonplane))
lab = {ports[s]: s for s in range(v)}
print('dim-1 entrances that are a source label', sum(1 for r in nonplane if dims[start[r]] == 1 and B[start[r]][0] in lab))
# who writes / reads them, per phase
wr = defaultdict(Counter); rd = defaultdict(Counter); rdcls = Counter(); wrcls = Counter()
for ph in order:
    if ph == 'centres': continue
    for e in c.get(ph) or []:
        if e[0] != 'add': continue
        a, b = e[1], e[2]
        if a in nonplane: wr[a][ph] += 1; wrcls[cls(b)] += 1
        if b in nonplane: rd[b][ph] += 1; rdcls[cls(a)] += 1
print('phases writing them', sum((wr[r] for r in nonplane), Counter()))
print('phases reading them', sum((rd[r] for r in nonplane), Counter()))
print('classes they read from', wrcls, ' classes that read them', rdcls)
# value replay in the hub's drop-free demand order (format-limits d1 demand531: precondition, A, B, cleanup; centres'
# decoder copies left out). Checked: targets end with exactly 5,824 nonzero source coefficients (= the hub's demand_nnz)
# and no entrance symbol (zero target response to entering content).
VORDER = ['precondition', 'A', 'B', 'cleanup']
val = [dict() for _ in range(n)]
for s in range(v): val[s] = {('x', s): Fr(1)}
for r in ent: val[r] = {('e', r): Fr(1)}
tread = Counter(); hread = Counter()
for ph in VORDER:
    for e in c.get(ph) or []:
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            if b in nonplane:
                kinds = tuple(sorted(Counter(key[0] for key in val[b]).items()))
                (tread if cls(a) == 'y' else hread)[(ph, kinds)] += 1
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
        elif e[0] == 'scale':
            _, a, num, den = e[:4]; k = Fr(num, den); val[a] = {key: x*k for key, x in val[a].items()}
print('targets with entrance symbols at the end:', sum(1 for t in range(v, 2*v) if any(k[0] == 'e' for k in val[t])),
      ' target nnz:', sum(sum(1 for k in val[t] if k[0] == 'x') for t in range(v, 2*v)))
print('reads of non-plane helpers by TARGETS (phase, (symbol kind, count) held):', dict(tread))
print('reads of non-plane helpers by others:', dict(hread))
