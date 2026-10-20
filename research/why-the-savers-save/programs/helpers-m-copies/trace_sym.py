import sys
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
syms = [int(x) for x in sys.argv[2].split(',')]
prog = program(w)
rows = [{q: F(1)} for q in range(n)]
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add':
        a, b = e[1], e[2]
        hit = [s for s in syms if s in rows[b] and s != a]
        if hit and not (a in syms):
            fd = None if e[-1] is None else len(basis(e[-1]))
            print(t, ph, i, 'reg', a, '+=', e[3], '/', e[4], '*', b, 'fdim', fd, {s: str(rows[b][s]) for s in hit}, '-> now', {s: str(rows[a].get(s, 0) + F(e[3], e[4])*rows[b][s]) for s in hit})
    apply(e, rows)
print('final rows containing syms:', {q: {s: str(rows[q][s]) for s in syms if s in rows[q]} for q in range(n) if any(s in rows[q] for s in syms) and q not in syms})
