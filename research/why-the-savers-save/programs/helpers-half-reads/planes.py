"""d-halfread: per packet plane of the U8 FULL-WORD: visitors, the gates at P, merge/copy roles, the readers of each
visitor after P (with the clean source content they take). Usage: python3 -I planes.py FULL-WORD.json.gz [out.json]"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); chain, bad = frame_replay(w, prog); assert not bad, bad[:3]
cls = lambda q: 'x' if q < v else ('y' if q < 2*v else 'h')
# values (all columns) and per-gate read values
rows = [{k: F(1)} for k in range(n)]; readval = {}
for T, (ph, e) in enumerate(prog):
    if e[0] == 'add': readval[T] = dict(rows[e[2]])
    apply(e, rows)
planes = defaultdict(list)
for q in range(n):
    for T, U in chain[q]:
        if len(U) == 3: planes[U].append((q, T))
print('planes', len(planes))
pats = Counter(); info = {}
for P, vis in sorted(planes.items()):
    gP = [(T, ph, e) for T, (ph, e) in enumerate(prog) if e[0] == 'add' and ph in ('A', 'B', 'tail') and basis(e[-1]) == P]
    info[P] = dict(vis=vis, gates=[(T, ph, e[1], e[2], e[3], e[4]) for T, ph, e in gP])
    pats[tuple(sorted(cls(q) for q, _ in vis)) + (len(gP),)] += 1
print(pats)
for P in sorted(planes)[:3]:
    print('plane', P, 'visitors', [(q, cls(q), T, prog[T][0] if T >= 0 else '-') for q, T in planes[P]])
    for T, ph, a, b, num, den in info[P]['gates']: print('   gate', T, ph, a, '+=', F(num, den), '*', b, ' b clean:', {k: str(x) for k, x in readval[T].items() if k < v})
# ---- roles and readers
def fmt(d): return {('x%d' % k if k < v else 'y%d' % (k-v) if k < 2*v else 'h%d' % k): str(x) for k, x in d.items()}
roles = {}
for P in sorted(planes):
    (Tm, phm, a, b, nm, dm), (Tc, phc, h3, a2, nc, dc) = info[P]['gates']
    assert a2 == a, P
    roles[P] = dict(merge=(Tm, phm, a, b, str(F(nm, dm))), copy=(Tc, phc, h3, str(F(nc, dc))))
def events_of(q):
    return [(T, ph, e) for T, (ph, e) in enumerate(prog) if q in regs(e)]
shown = 0
for P in sorted(planes):
    r = roles[P]; h3 = r['copy'][2]
    if shown >= 4: break
    if P not in [sorted(planes)[i] for i in (0, 1, 2, 5)]: continue
    shown += 1
    print('\nPLANE', P, r)
    for T, ph, e in events_of(h3):
        if e[0] != 'add': print('  ', T, ph, e); continue
        rd = 'READ by %s%d' % (cls(e[1]), e[1]) if e[2] == h3 else 'WRITE from %s%d' % (cls(e[2]), e[2])
        print('  ', T, ph, rd, 'q', F(e[3], e[4]), 'frame dim', len(basis(e[-1])) if e[-1] is not None else None,
              ' h3 val' if e[2] == h3 else ' src val', fmt(readval[T]) if len(readval[T]) < 14 else '(%d terms)' % len(readval[T]))
