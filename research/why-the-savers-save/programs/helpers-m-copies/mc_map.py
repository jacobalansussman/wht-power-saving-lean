"""d-mcopy-code STEP 1: exact map of the 112 M copies (entering helpers that never visit a rank-3 frame) of the U8 word.
Usage: python3 -I mc_map.py FULL-WORD.json.gz mc_map.json"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
start, end, chain, paid, runs = frames(w)
paid[8] += 8
print('paid (with 8 centre children):', dict(sorted(paid.items())), '== report' , dict(sorted(paid.items())) == {int(k): a for k, a in w['report']['local_profile'].items()})
vis3 = {q for q in range(n) if any(len(U) == 3 for _, U in chain[q])}
ent = {int(q) for q in w['start_extra']}
M = sorted(q for q in range(2*v, n) if q in ent and q not in vis3)
print('M copies', len(M), 'entrance dims', Counter(len(start[q]) for q in M), 'widths', Counter(len(end[q]) - len(start[q]) for q in M))
cls = lambda q: 'x' if q < v else ('y' if q < 2*v else ('M' if q in M else 'h'))
prog = program(w)
# symbolic replay: every register starts as its own symbol
rows = [{q: F(1)} for q in range(n)]
ev = defaultdict(list)   # q -> events
def summ(d):
    c = Counter(cls(k) if k != None else '?' for k in d)
    return dict(c)
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add':
        a, b = e[1], e[2]; k = F(e[3], e[4]); U = e[-1]
        fd = None if U is None else len(basis(U))
        if a in M or b in M:
            val_b = dict(rows[b])
            for q, role, p in ((a, 'w', b), (b, 'r', a)):
                if q in M: ev[q].append(dict(t=t, ph=ph, i=i, role=role, partner=p, pcls=cls(p), k=str(k), fdim=fd, U=list(basis(U)) if U is not None else None,
                                              val_b={str(x): str(y) for x, y in val_b.items()}))
        apply(e, rows)
    else:
        if e[1] in M: ev[e[1]].append(dict(t=t, ph=ph, i=i, role='s', k=str(F(e[2], e[3]))))
        apply(e, rows)
# endpoint check
bad = [q for q in range(n) if (rows[q] != ({q: F(1), q - v: F(1)} if v <= q < 2*v and False else rows[q]))]
# summary per M copy
out = {}
for q in M:
    E = ev[q]
    reads = [x for x in E if x['role'] == 'r']; writes = [x for x in E if x['role'] == 'w']
    out[q] = dict(entrance=list(start[q]), end=list(end[q]), chain=[[str(tag), list(U)] for tag, U in chain[q]],
                  nreads=Counter(x['pcls'] + ':' + x['ph'] for x in reads), nwrites=Counter(x['pcls'] + ':' + x['ph'] for x in writes),
                  events=E)
json.dump({str(q): o for q, o in out.items()}, open(sys.argv[2], 'w'))
# print 3 examples compactly, one of each entrance dim
shown = set()
for q in M:
    d = len(start[q])
    if d in shown: continue
    shown.add(d)
    print('\n=== M copy', q, 'entrance', start[q], 'dim', d, 'end dim', len(end[q]))
    print('  chain dims', [len(U) for _, U in chain[q]], [str(tag) for tag, U in chain[q]])
    for x in ev[q]:
        if x['role'] == 's': print('  ', x['t'], x['ph'], 'scale', x['k']); continue
        vb = x['val_b']
        sym = Counter(cls(int(k)) for k in vb)
        print('  ', x['t'], x['ph'], x['i'], x['role'], x['pcls'], x['partner'], 'k', x['k'], 'fdim', x['fdim'], 'read value classes', dict(sym) if x['role']=='r' else '', 'partner-val' if x['role']=='w' else '', dict(sym) if x['role']=='w' else '')
