#!/usr/bin/env python3
"""d-bit-exits / totals_exit.py: can #589's nine retained totals (never read by a z register) be restored right after
the scatter at their own hyperplane q-perp (an exit at dim 7), before B moves their builders away?
Tests at tau = just after the scatter, for each total h at U = its frame: (i) no register outside U holds e_h;
(ii) e_h lies in the F2 span of the values of registers whose frame is inside U, and lists how many of a minimal-ish
solution's registers sit strictly below U (they would have to climb to U, splitting a run).
usage: python3 -I totals_exit.py WORD.json.gz"""
import gzip, json, os, sys
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames_lib import Frames
w = json.load(gzip.open(sys.argv[1], 'rt')); V, R = w['v'], w['R']; N = 2*V+R; F = Frames(w)
cur = list(w['start']); val = [1 << i for i in range(N)]
for sec in ('D', 'A'):
    for f, t, ss in w[sec]:
        for r in [t, *ss]: cur[r] = f
        for s in ss: val[t] ^= val[s]
snap = [val[q] for q in w['ret']]
for t, ks in enumerate(w['scatter']):
    for k in ks: val[V+t] ^= snap[k]
res = []
for h in w['ret']:
    U = cur[h]
    outside = [r for r in range(N) if val[r] >> h & 1 and not F.le(cur[r], U)]
    inside = [r for r in range(N) if F.le(cur[r], U) and r != h and not (V <= r < 2*V)]
    target = val[h] ^ (1 << h)          # what h must lose
    basis = {}
    for r in inside:
        x = val[r]; s = {r}
        while x:
            b = x.bit_length()-1
            if b in basis: x ^= basis[b][0]; s ^= basis[b][1]
            else: basis[b] = (x, s); break
    x = target; sol = set(); ok = True
    while x:
        b = x.bit_length()-1
        if b not in basis: ok = False; break
        x ^= basis[b][0]; sol ^= basis[b][1]
    below = [r for r in sol if cur[r] != U] if ok else None
    res.append(dict(total=h, U_dim=F.dim[U], dirt_outside=len(outside), span_ok=ok, solution_size=len(sol) if ok else None,
                    solution_below_U=len(below) if ok else None))
print(json.dumps(res, indent=1))
