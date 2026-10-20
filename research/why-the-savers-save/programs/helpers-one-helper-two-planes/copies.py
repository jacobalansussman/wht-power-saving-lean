#!/usr/bin/env python3
"""Step 1b: the copy-holder roles. For each plane, the 2 visitors leaving P with the packet (copy holders) and the
proper-row visitor. For each copy holder: all gates where it is READ (others += c*q) between leaving P and the start of
cleanup, with reader class and gate frame dims; and all gates where it is WRITTEN after P (content changes).
Usage: python3 -I copies.py FULL-WORD.json.gz survey.json copies.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
from collections import Counter, defaultdict
w = load(sys.argv[1]); S = json.load(open(sys.argv[2])); v, h, R, n, full, start, end = geometry(w)
segs, comp = program(w); off = S['seg_offsets']
gl = []
for name, g in segs: gl += [(name, e) for e in g]
def cls(q):
    if q < v: return 'x'
    if q < 2 * v: return 'y'
    hh = q - 2 * v
    return 'fac' if hh < 224 else ('face' if hh < 392 else 'M')
byreg = defaultdict(list)
for t, (name, e) in enumerate(gl):
    for q in registers(e): byreg[q].append(t)
out = []; rc = Counter(); nreads = Counter(); pr = Counter()
for r in S['planes']:
    P = tuple(r['plane'])
    for d in r['visitors']:
        q = d['q']; role = 'copy' if d['out_spans_P'] else 'proper'
        reads = []; writes = []
        for t in byreg[q]:
            if t <= d['t_in']: continue
            name, e = gl[t]
            if name in ('cleanup', 'compensation', 'precondition'): continue
            fr = e[-1]
            if e[0] == 'add' and e[2] == q: reads.append((t, name, e[1], cls(e[1]), len(basis(fr)) if isinstance(fr, list) else None))
            elif e[1] == q: writes.append((t, name, e[2] if e[0] == 'add' else None, cls(e[2]) if e[0] == 'add' else 'scale', len(basis(fr)) if isinstance(fr, list) else None))
        out.append(dict(plane=list(P), block=r['block'], q=q, cls=d['cls'], role=role, t_in=d['t_in'], t_out=d['t_out'], reads=reads, writes=writes))
        if role == 'copy':
            rc[(d['cls'], tuple(sorted(Counter(x[3] for x in reads).items())))] += 1
            nreads[len(reads)] += 1
        else: pr[d['cls']] += 1
json.dump(out, open(sys.argv[3], 'w'))
print('copy holders by (class, reader classes before cleanup):')
for k, x in rc.most_common(): print('  ', k, x)
print('copy holders by number of non-cleanup reads', dict(nreads)); print('proper-row visitors by class', dict(pr))
