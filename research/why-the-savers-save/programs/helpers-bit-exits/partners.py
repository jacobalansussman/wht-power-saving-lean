#!/usr/bin/env python3
"""d-bit-exits / partners.py: for every #589 helper, can its cleanup partners all sit below one hyperplane U that
contains its end-of-B frame? This is a necessary condition for a hyperplane exit that keeps the literal cleanup gates
(in any order): a C gate runs at a common frame containing every operand's current frame, and frames only grow.
usage: python3 -I partners.py WORD.json.gz [out.json]"""
import gzip, json, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames_lib import Frames

w = json.load(gzip.open(sys.argv[1], 'rt'))
V, R = w['v'], w['R']; N = 2 * V + R
F = Frames(w)
cur = list(w['start'])
for sec in ('D', 'A', 'B'):
    for f, t, ss in w[sec]:
        for r in [t, *ss]: cur[r] = f
endB = list(cur)
partners = defaultdict(set)
for f, t, ss in w['C']:
    regs = [t, *ss]
    for r in regs:
        for q in regs:
            if q != r: partners[r].add(q)
kind = lambda r: 'x' if r < V else ('y' if r < 2 * V else 'h')
res = Counter(); detail = Counter(); outside = Counter()
for h in range(2 * V, N):
    Fh = endB[h]; dh = F.dim[Fh]
    j = F.join([endB[q] for q in partners[h]] + [Fh])
    res['%d->%d' % (dh, F.dim[j])] += 1
    bad = [q for q in partners[h] if not F.le(endB[q], Fh)]
    detail['end dim %d: partners outside own frame: %s' % (dh, '+'.join(sorted({kind(q) for q in bad})) or 'none')] += 1
    outside[len(bad)] += 1
out = dict(join_dims=dict(sorted(res.items())), outside_kinds=dict(detail),
           outside_partner_count=dict(sorted(outside.items())))
print(json.dumps(out, indent=1))
if len(sys.argv) > 2: json.dump(out, open(sys.argv[2], 'w'), indent=1)
