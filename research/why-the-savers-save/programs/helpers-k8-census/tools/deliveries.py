#!/usr/bin/env python3
"""k8-pairs: how each pm.py need (target t, block b with |demand_t & b| = 4) is delivered: the multiset of
|content & b| over t's reads (precondition/B/tail; clean replay over sources, cut/compensation skipped, centres
ignored) whose content meets b, with the reading register's class.  Need kind: 'own' if b is t's own block, else 'face'.
Usage: python3 -I deliveries.py WORD U8_BLOCKS DEMAND out.json"""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
blocks = json.load(open(sys.argv[2]))['blocks']; blk = {s: i for i, x in enumerate(blocks) for s in x}
d = json.load(open(sys.argv[3])); dem = {int(t): {int(s) for s in dd} for t, dd in d['demand'].items()}
v = c['v']; R = c['R']; n = 2*v + R
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: Fr(1)}
rd = defaultdict(list)   # (t, b) -> list of (cls, |content & b|, |content|)
for ph in ['A', 'precondition', 'B', 'tail', 'cleanup']:
    for e in c.get(ph) or []:
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            if v <= a < 2*v and val[b]:
                fp = Counter(blk[s] for s in val[b])
                for bb, x in fp.items(): rd[(a - v, bb)].append((ph[0] + cls(b), x, len(val[b])))
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
        else:
            _, a, num, den = e[:4]; k = Fr(num, den); val[a] = {key: x*k for key, x in val[a].items()}
per = defaultdict(Counter)
for t in range(v):
    for b in range(len(blocks)):
        S = dem[t] & set(blocks[b])
        if len(S) != 4: continue
        kind = 'own' if blk[t] == b else 'face'
        sig = tuple(sorted((cl, x, tot if tot <= 4 else '>4') for cl, x, tot in rd[(t, b)]))
        per[(b, kind)][sig] += 1
out = {}
for (b, kind), cn in sorted(per.items()):
    out[f'{b}:{kind}'] = {str(k): x for k, x in cn.most_common()}
json.dump(out, open(sys.argv[4], 'w'), indent=0)
for (b, kind), cn in sorted(per.items()):
    if b in (0, 1, 7):
        print(b, kind, sum(cn.values()), cn.most_common(6))
