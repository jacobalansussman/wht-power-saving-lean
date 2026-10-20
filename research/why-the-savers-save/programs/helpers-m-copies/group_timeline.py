"""print one group's complete timeline: gates touching its M copies, and gates among/into its 4 targets in B, with values"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
G = json.load(open('mc_groups.json'))[int(sys.argv[2])]
mem = set(G['members']); T = set(G['targets'])
start, end, chain, paid, runs = frames(w)
ports = w['ports']
prog = program(w)
rows = [{q: F(1)} for q in range(n)]
def show(d):
    xs = sorted(s for s in d if s < v); es = sorted(s for s in d if s in mem); hs = sorted(s for s in d if s >= 2*v and s not in mem); ys = sorted(s for s in d if v <= s < 2*v)
    return f"x{len(xs)} M{[ (s, str(d[s])) for s in es]} h{len(hs)} y{ys}"
print('members', sorted(mem), 'entrances', {q: start[q] for q in mem})
print('targets', sorted(T), 'ports', {t: ports[t - v] for t in T}, 'end frames u^perp dims', {t: len(end[t]) for t in T})
for t, (ph, i, e) in enumerate(prog):
    if e[0] != 'add': 
        if e[1] in mem or e[1] in T: 
            if ph in ('B','A','tail'): print(t, ph, 'scale', e[1], e[2], e[3])
        continue
    a, b = e[1], e[2]
    rel = (a in mem or b in mem) or (ph in ('B', 'A', 'tail') and a in T and (b in T or b in mem))
    if rel and ph != 'cleanup':
        fd = None if e[-1] is None else len(basis(e[-1]))
        print(t, ph, i, f"{a}{'*' if a in mem else ('T' if a in T else '')} += {F(e[3], e[4])} * {b}{'*' if b in mem else ('T' if b in T else '')}", 'frame', fd, list(basis(e[-1])) if e[-1] is not None else None, '| val', show(rows[b]))
    apply(e, rows)
