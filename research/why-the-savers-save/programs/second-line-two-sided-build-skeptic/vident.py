#!/usr/bin/env python3
"""v-an-build: the invocation identity of the U8X word, own code.
Forward word F = g8(Mi) . M . g1(M) (Lean InvThm: g1 at the start frames, main phase, g8 at the end), exactly in Fractions,
must be (x, y, s) -> (x, y + K x, s) with K = M[y, x]; then the same word run on raw complex arrays (f = 1..8 columns,
random Gaussian-integer data) against the formula.  usage: python3 -I vident.py WORD.json.gz"""
import sys, random, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vreplay import load, forward, inverse, gate_parts, Rows
from fractions import Fraction as F
d = load(sys.argv[1]); v, R = d['v'], d['R']; n = 2 * v + R
M = forward(d); Mi = inverse(d)
def run(state, addf, negf):
    """apply the forward word to a state given by callbacks (state[r] = value of register r)."""
    pass
# exact composite on rows
X = Rows(n)
for T in range(v):                                       # g1: y_T -= sum_q M[y_T, q] s_q   (slots q)
    for q, c in M[v + T].items():
        if q >= 2 * v: X.add(v + T, q, -c)
for ph in ('A',):
    for g in d[ph]:
        _, adds, ng = gate_parts(g)
        for t, s, c in adds: X.add(t, s, c)
copies = {k: dict(X.row[src]) for k, src, _ in d['ret']}
for T, ent in enumerate(d['scat']['table']):
    for k, a, b in ent: X.addform(v + T, copies[k], F(a, b))
for g in d['B']:
    _, adds, ng = gate_parts(g)
    for t, s, c in adds: X.add(t, s, c)
new = {}
for q in range(2 * v, n):                                # g8: slot q := sum_k Mi[q, k] z_k
    acc = {}
    for k, c in Mi[q].items():
        for j, val in X.row[k].items():
            acc[j] = acc.get(j, 0) + c * val
    new[q] = {j: x for j, x in acc.items() if x}
for q, r in new.items(): X.row[q] = r
okx = all(X.row[r] == {r: F(1)} for r in range(v))
oks = all(X.row[q] == {q: F(1)} for q in range(2 * v, n))
oky = all(X.row[v + T].get(v + T) == 1 and all((k < v) or k == v + T for k in X.row[v + T]) for T in range(v))
K = [{k: c for k, c in M[v + T].items() if k < v} for T in range(v)]
okK = all({k: c for k, c in X.row[v + T].items() if k < v} == K[T] for T in range(v))
print('exact composite: x fixed', okx, '| slots restored', oks, '| y -> y + (x terms only)', oky, '| x-block = K = M[y,x]', okK)
ents = sorted(set(c for T in range(v) for c in K[T].values()))
print('K: nonzeros per row', sorted(set(len(K[T]) for T in range(v))), 'entry values', [str(e) for e in ents])
# raw complex arrays
random.seed(15)
allok = True
for f in range(1, 9):
    z = [[complex(random.randint(-9, 9), random.randint(-9, 9)) for _ in range(f)] for _ in range(n)]
    z0 = [list(col) for col in z]
    def addv(t, s, c):
        cc = complex(c)
        z[t] = [a + cc * b for a, b in zip(z[t], z[s])]
    Mys = [[(q, complex(c)) for q, c in M[v + T].items() if q >= 2 * v] for T in range(v)]
    for T in range(v):
        for q, c in Mys[T]: z[v + T] = [a - c * b for a, b in zip(z[v + T], z[q])]
    for g in d['A']:
        for t, s, c in gate_parts(g)[1]: addv(t, s, c)
    cp = {k: list(z[src]) for k, src, _ in d['ret']}
    for T, ent in enumerate(d['scat']['table']):
        for k, a, b in ent:
            c = a / b; z[v + T] = [x + c * y for x, y in zip(z[v + T], cp[k])]
    for g in d['B']:
        for t, s, c in gate_parts(g)[1]: addv(t, s, c)
    snap = [list(col) for col in z]
    for q in range(2 * v, n):
        acc = [0j] * f
        for k, c in Mi[q].items():
            cc = complex(c); acc = [a + cc * b for a, b in zip(acc, snap[k])]
        z[q] = acc
    err = 0.0
    for r in range(v): err = max(err, max(abs(a - b) for a, b in zip(z[r], z0[r])))
    for q in range(2 * v, n): err = max(err, max(abs(a - b) for a, b in zip(z[q], z0[q])))
    for T in range(v):
        want = list(z0[v + T])
        for k, c in K[T].items():
            cc = complex(c); want = [a + cc * b for a, b in zip(want, z0[k])]
        err = max(err, max(abs(a - b) for a, b in zip(z[v + T], want)))
    print('raw complex arrays f = %d: max error %.3e' % (f, err))
    allok &= err < 1e-6
print('ALL raw checks pass:', allok)
