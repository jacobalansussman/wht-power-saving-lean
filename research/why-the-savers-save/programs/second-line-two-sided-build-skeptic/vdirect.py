#!/usr/bin/env python3
"""v-an-build: count gates whose register set holds both an x (r < v) and a y (v <= r < 2v) role, per phase.
Own code; gate fields as in gx.py gate_regs.  usage: python3 -I vdirect.py WORD.json.gz ..."""
import gzip, json, sys
from collections import Counter
def regs(g):
    if g[0] == 'out': return [g[2]] + [t for t, _, _ in g[3]] + list(g[4])
    if g[0] == 'in': return [g[2]] + [s for s, _, _ in g[3]]
    if g[0] == 'neg': return [g[2]]
    raise ValueError(g[0])
for p in sys.argv[1:]:
    d = json.load(gzip.open(p, 'rt')); v = d['v']
    res = {}
    for ph in ('A', 'B'):
        n = 0; kinds = Counter()
        for g in d[ph]:
            rs = regs(g)
            if any(r < v for r in rs) and any(v <= r < 2 * v for r in rs):
                n += 1; kinds[g[0]] += 1
        res[ph] = (n, dict(kinds))
    print(p.split('/')[-1], 'format', d['format'], 'v', v, 'R', d['R'], 'gates A/B', len(d['A']), len(d['B']), 'direct x-y gates', res)
