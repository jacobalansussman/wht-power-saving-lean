#!/usr/bin/env python3
"""Print the full life of given registers: every gate touching them (time, phase, role, partner, gate frame, frame
climb, content after), in the audit's program order.  Usage: python3 -I life.py FULL-WORD.json.gz q1 q2 ..."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
w = load(sys.argv[1]); v, h, R, n, full, start, end = geometry(w); ports = w['ports']
Q = [int(a) for a in sys.argv[2:]]
segs, comp = program(w)
rows = [{k: F(1)} for k in range(n)]
def fmt(r):
    return ' '.join(('x%d' % k if k < v else 'y%d' % k if k < 2*v else 'd%d' % k) + ('' if x == 1 else '*' + str(x)) for k, x in sorted(r.items()))
cur = list(start); t = 0
for q in Q: print('reg', q, 'start', start[q], 'end', end[q], 'label' if q < v else '', ports[q] if q < v else '')
for name, g in segs:
    for e in g:
        fr = e[-1] if (name not in ('precondition', 'compensation', 'cut', 'cleanup')) else None
        before = {q: cur[q] for q in registers(e)}
        if fr is not None and name in ('A', 'B', 'tail'):
            for q in registers(e): cur[q] = basis(fr)
        apply(e, rows)
        hit = [q for q in registers(e) if q in Q]
        if hit:
            role = 'W' if e[1] in Q else 'R'
            other = [q for q in registers(e) if q not in Q] or registers(e)
            print('%5d %-12s %s %s' % (t, name, ' '.join('%s%d' % ('W' if q == e[1] else 'r', q) for q in registers(e)), e[3:5] if e[0]=='add' else e[2:4]),
                  'F', list(basis(e[-1])) if isinstance(e[-1], list) else e[-1], '|', ' ; '.join('%d:%s' % (q, fmt(rows[q])) for q in registers(e) if q in Q))
        t += 1
