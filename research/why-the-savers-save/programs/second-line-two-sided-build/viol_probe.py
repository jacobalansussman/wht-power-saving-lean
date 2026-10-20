#!/usr/bin/env python3
"""an-build: where the pass-1 violations of two-sided slack sit (own code). usage: python3 -I viol_probe.py FILE"""
import gzip, json, sys
from collections import Counter
d = json.load(gzip.open(sys.argv[1]))
h, v, R = d['h'], d['v'], d['R']; n = 2 * v + R
A, B = d['A'], d['B']; GA = A + B; nA = len(A)
tot = [s_ for _, s_, _ in d['ret']]
def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
RG = [regs(g) for g in GA]
def cl(r): return 'x' if r < v else 'y' if r < 2 * v else ('t' if r in tot else 's')
s = [False] * n
for r in range(v, 2 * v): s[r] = True
G = []; first_s = {}
for i, rs in enumerate(RG):
    if i == nA:
        for t in tot: s[t] = True; first_s.setdefault(t, i)
    j = any(s[r] for r in rs)
    if j:
        for r in rs:
            if not s[r]: first_s[r] = i
            s[r] = True
    G.append(j)
a = [False] * n
for r in range(v): a[r] = True
Aat = [False] * len(RG)
for i in range(len(RG) - 1, -1, -1):
    j = any(a[r] for r in RG[i])
    if j:
        for r in RG[i]: a[r] = True
    Aat[i] = j
bad = [i for i in range(len(RG)) if G[i] and Aat[i]]
print(sys.argv[1].split('/')[-1], 'violating gates', len(bad), '(A:', sum(i < nA for i in bad), 'B:', sum(i >= nA for i in bad), ')')
print('  kinds of violating gates (classes of registers named):', Counter(''.join(sorted(cl(r) for r in RG[i])) for i in bad).most_common(8))
print('  x registers named in a violating gate:', len({r for i in bad for r in RG[i] if r < v}))
# do x's meet totals after the cut?  do x's meet y's?
print('  B gates naming an x and a total:', sum(1 for g in RG[nA:] if any(r < v for r in g) and any(r in tot for r in g)),
      '| naming an x and a y:', sum(1 for g in RG[nA:] if any(r < v for r in g) and any(v <= r < 2 * v for r in g)))
# first violating gate: how did s reach it?
i0 = bad[0] if bad else None
if i0 is not None:
    print('  first violating gate', i0, GA[i0][:3], 'regs', [(r, cl(r), first_s.get(r)) for r in RG[i0]][:8])
