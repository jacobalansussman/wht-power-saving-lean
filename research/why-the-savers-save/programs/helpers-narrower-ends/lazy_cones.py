#!/usr/bin/env python3
"""d-hpe-multigate STEP 1c: how far can end-B frames fall if only frame annotations change (scalar gates and order fixed)?
Replays A and B with every gate at the JOIN of its operands' current frames (the least legal frame; centres still forced
to their interface frames before B, targets at 0 before B), then recomputes the two cones of cones.py with these least
end-B frames.  Umin_lazy(q) is a lower bound on the end of helper q for ANY frame re-optimisation of A/B and ANY cleanup
with the same prefix scalars (the cones are frame-free; frames only enter through endB).
Usage: python3 -I lazy_cones.py FULL-WORD.json.gz cones.json out.json"""
import sys, gzip, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from cones_lib import basis, contained, registers, apply, invert
w = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read())); cones = json.load(open(sys.argv[2]))['helpers']
v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v+R
full = basis([1 << i for i in range(h)])
start = [(p,) for p in w['ports']]+[()]*(v+R)
for q, U in w['start_extra'].items(): start[int(q)] = basis(U)
cur = start.copy(); lazyB = 0; tighter = Counter()
for section in ('A', 'B'):
    if section == 'B':
        for q, U in zip(w['centers'], w['center_frames']):
            U = basis(U); assert contained(cur[q], U); cur[q] = U
    for e in w[section]:
        rs = registers(e)
        J = basis([x for q in rs for x in cur[q]])
        assert contained(J, basis(e[-1]))
        if J != basis(e[-1]): tighter[(section, len(basis(e[-1])) - len(J))] += 1
        for q in rs: cur[q] = J
endB_lazy = cur.copy()
# the word's own end-B frames
cur = start.copy()
for section in ('A', 'B'):
    for e in w[section]:
        for q in registers(e): cur[q] = basis(e[-1])
endB = cur.copy()
print('gates whose frame is above the join of their operands (section, excess dim):', dict(sorted(tighter.items())))
print('registers whose end-B frame falls under lazy frames:', Counter((('x' if q < v else 'y' if q < 2*v else 'h'), len(endB[q]), len(endB_lazy[q])) for q in range(n) if endB[q] != endB_lazy[q]))
ident = lambda: [{k: F(1)} for k in range(n)]
main = w['A']+w['cut']+w['B']
raw = ident()
for e in main: apply(e, raw)
dirt = [{} for _ in range(n)]
for y in range(v, 2*v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2*v}
for e in w['precondition']: apply(e, dirt)
compensation = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2*v) for k, a in sorted(dirt[y].items())]
prefix = w['precondition']+compensation+main
p = ident(); pinv = ident()
for e in prefix: apply(e, p)
for e in reversed(prefix): apply(invert(e), pinv)
col = defaultdict(set)
for k in range(n):
    for q in p[k]: col[q].add(k)
out = {}; summ = Counter()
for q in range(2*v, n):
    rowk = set(pinv[q]) | {q}; colk = col[q] | {q}
    Ur = basis([x for k in rowk for x in endB_lazy[k]]); Uc = basis([x for k in colk for x in endB_lazy[k]])
    U = basis(list(Ur)+list(Uc))
    out[q] = dict(end_dim=cones[str(q)]['end_dim'], Umin_dim=cones[str(q)]['Umin_dim'], Ulazy_dim=len(U), Ulazy_row=len(Ur), Ulazy_col=len(Uc), Ulazy=list(U))
    summ[(out[q]['end_dim'], len(U))] += 1
print('(current end dim, least end dim under lazy A/B frames) : helpers', dict(sorted(summ.items())))
json.dump(dict(summary={f'{a}->{b}': x for (a, b), x in sorted(summ.items())}, helpers=out,
               endB_lazy=[list(f) for f in endB_lazy]), open(sys.argv[3], 'w'), indent=0)
