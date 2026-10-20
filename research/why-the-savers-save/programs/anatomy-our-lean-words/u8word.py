#!/usr/bin/env python3
"""an-ours: class-by-class climbs of the U8 local word as lowered by us (k-lower/out/U8X-parsed.json.gz, gcert/1, our
own data file; this U8X version carries the 112 source-copy helpers, R = 616).  Own replay, after the semantics of our
checker tools/gx/gxcore.py chk_replay (x = sources 0..v-1, y = targets v..2v-1, s = slots, c = retained totals)."""
import gzip, json, math, sys
from collections import Counter
d = json.load(gzip.open(sys.argv[1]))
h, v, R = d['h'], d['v'], d['R']; fr = d['frames']; dim = [len(f) for f in fr]
n = 2 * v + R
cls = lambda r: 'x' if r < v else 'y' if r < 2 * v else 's'
cur = list(d['start']); H = {k: Counter() for k in 'xysc'}
def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
def climb(r, f):
    if cur[r] != f:
        assert dim[f] > dim[cur[r]]
        H[cls(r)][dim[f] - dim[cur[r]]] += 1; cur[r] = f
for i, g in enumerate(d['A'] + d['B']):
    if i == len(d['A']):
        for _, s, f in d['ret']:
            assert cur[s] == f; H['c'][dim[f]] += 1
    for r in regs(g): climb(r, g[1])
for r in range(n): climb(r, d['final'][r])
ext = dict((r, f) for r, f in d['ext'])
ent = sum(dim[f] for f in ext.values())
rank = {k: sum(a * b for a, b in w.items()) for k, w in H.items()}
nm = {k: sum(a * b * math.log(73 / a) for a, b in w.items()) for k, w in H.items()}
tot = sum(rank.values())
print('h', h, 'v', v, 'R', R, 'registers', n)
for k in 'xysc':
    print(f"  class {k}: blocks {dict(sorted(H[k].items()))}  rank {rank[k]}  narrow(73) {nm[k]:.1f}  ({100*nm[k]/sum(nm.values()):.1f} % of narrow)")
print('x starts on its port line:', all(fr[d['start'][t]] == [d['ports'][t]] for t in range(v)),
      '| y ends on a hyperplane:', all(dim[d['final'][v + t]] == h - 1 for t in range(v)))
print('total rank', tot, '= R h + 2v(h-1) + cst - ent =', R * h + 2 * v * (h - 1) + rank['c'] - ent, '(cst', rank['c'], 'ent', ent, ')')
allH = Counter()
for w in H.values(): allH.update(w)
print('whole local profile', dict(sorted(allH.items())))
json.dump(dict(H={k: dict(w) for k, w in H.items()}, rank=rank, narrow=nm, ent=ent), open('u8word.json', 'w'), indent=1)
