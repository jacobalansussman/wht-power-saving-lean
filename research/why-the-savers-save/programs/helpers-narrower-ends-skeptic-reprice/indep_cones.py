"""Skeptic's independent cone census (own code; mod-prime arithmetic; product check of the inverse)."""
import sys, gzip, json
from collections import Counter
PR = (1 << 61) - 1
def md(num, den): return num % PR * pow(den % PR, PR - 2, PR) % PR
def rank2(vs):
    piv = {}
    for a in vs:
        while a:
            b = a.bit_length() - 1
            if b in piv: a ^= piv[b]
            else: piv[b] = a; break
    return piv
def span(vs): return frozenset(rank2(vs).values())
def dim(vs): return len(rank2(vs))
def regs(e): return e[1:3] if e[0] == 'add' else e[1:2]
def app(e, rows):  # left row operation, mod PR
    if e[0] == 'add':
        t, s = e[1], e[2]; q = md(e[3], e[4]); d = rows[t]
        for k, x in rows[s].items():
            y = (d.get(k, 0) + q * x) % PR
            if y: d[k] = y
            else: d.pop(k, None)
    else:
        t = e[1]; q = md(e[2], e[3]); rows[t] = {k: q * x % PR for k, x in rows[t].items()}
def inv(e):
    if e[0] == 'add': return ['add', e[1], e[2], -e[3], e[4]]
    return ['scale', e[1], e[3], e[2]]
w = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2 * v + R
I = lambda: [{k: 1} for k in range(n)]
main = w['A'] + w['cut'] + w['B']
# compensation per hub rule, exact rationals needed for the gate scalars -> use fractions here only
from fractions import Fraction as F
def appF(e, rows):
    if e[0] == 'add':
        t, s = e[1], e[2]; q = F(e[3], e[4]); d = rows[t]
        for k, x in rows[s].items():
            y = d.get(k, 0) + q * x
            if y: d[k] = y
            else: d.pop(k, None)
    else:
        t = e[1]; q = F(e[2], e[3]); rows[t] = {k: q * x for k, x in rows[t].items()}
raw = [{k: F(1)} for k in range(n)]
for e in main: appF(e, raw)
dirt = [{} for _ in range(n)]
for y in range(v, 2 * v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2 * v}
for e in w['precondition']: appF(e, dirt)
comp = [['add', y, k, -a.numerator, a.denominator] for y in range(v, 2 * v) for k, a in sorted(dirt[y].items())]
prefix = w['precondition'] + comp + main
P = I()
for e in prefix: app(e, P)
Q = I()
for e in reversed(prefix): app(inv(e), Q)
# product check Q*P == I  (row i of Q*P = sum_k Q[i][k] * P[k])
bad = 0
for i in range(n):
    acc = {}
    for k, a in Q[i].items():
        for j, b in P[k].items():
            acc[j] = (acc.get(j, 0) + a * b) % PR
    acc = {j: x for j, x in acc.items() if x}
    if acc != {i: 1}: bad += 1
print('Q*P == I rows failing:', bad)
start = [[p] for p in w['ports']] + [[] for _ in range(v + R)]
for q, U in w['start_extra'].items(): start[int(q)] = list(U)
def frames(lazy):
    cur = [span(s) for s in start]
    for sec in ('A', 'B'):
        if sec == 'B':
            for q, U in zip(w['centers'], w['center_frames']):
                assert cur[q] <= span(list(cur[q]) + U) and dim(list(cur[q]) + U) == dim(U); cur[q] = span(U)
        for e in w[sec]:
            rs = regs(e)
            J = span([x for q in rs for x in cur[q]])
            Ue = span(e[-1])
            assert dim(list(J) + list(Ue)) == len(Ue), 'gate frame below join'
            for q in rs: cur[q] = J if lazy else Ue
    return cur
colsup = [set() for _ in range(n)]
for k in range(n):
    for j in P[k]: colsup[j].add(k)
fin = {int(q): U for q, U in w['helper_final'].items()}
res = {}
for lazy in (False, True):
    eB = frames(lazy); tab = Counter(); rowonly = Counter()
    for q in range(2 * v, n):
        rowk = set(Q[q]) | {q}; colk = colsup[q] | {q}
        dr = dim([x for k in rowk for x in eB[k]]); dc = dim([x for k in colk for x in eB[k]])
        du = dim([x for k in rowk | colk for x in eB[k]])
        enddim = dim(fin.get(q, [1 << i for i in range(h)]))
        tab[(enddim, du)] += 1
        if enddim == 9: rowonly['row9' if dr == 9 else 'col-needed']+= 1
        res.setdefault(q, {})['lazy' if lazy else 'fixed'] = du
    print('lazy' if lazy else 'fixed endB', dict(sorted(tab.items())), dict(rowonly))
json.dump({str(q): r for q, r in res.items()}, open(sys.argv[2], 'w'))
