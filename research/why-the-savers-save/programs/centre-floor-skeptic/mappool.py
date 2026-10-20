#!/usr/bin/env python3
"""v-centre / mappool.py: map ctr-window's pool (its own line indexing) to my lines by root vectors (taken from
cx-medium's designs_cx.design('P112'), data only), check each pool set is a basis of R^8 (exact Fraction determinant
via my own elimination), that the first 1260 are orthogonal frames, and write pool_mine.txt (my indices) + integer
weights c_S = round(mu_S * 2^30)."""
import sys, os, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from designs_cx import design
theirs = design('P112')
mine = [tuple(v) for v in json.load(open(os.path.join(HERE, 'data.json')))['lines']]
def canon(v):
    v = list(v)
    # scale to integer vector proportional, sign: first nonzero positive, then reduce by gcd
    from math import gcd
    fr = [Fr(x) for x in v]
    den = 1
    for x in fr: den = den * x.denominator // gcd(den, x.denominator)
    iv = [int(x * den) for x in fr]
    g = 0
    for x in iv: g = gcd(g, abs(x))
    iv = [x // g for x in iv]
    for x in iv:
        if x: 
            if x < 0: iv = [-y for y in iv]
            break
    return tuple(iv)
# use THEIR roots directly; my own labels from them (F2 row space of the dense pairing, own elimination)
mine = [tuple(v) for v in theirs]
mp = list(range(112))
ipf = lambda a, b: sum(Fr(x) * Fr(y) for x, y in zip(a, b))
nrm = ipf(mine[0], mine[0])
rows = []
for t in range(112):
    r = 0
    for u in range(112):
        x = ipf(mine[t], mine[u]) / (nrm / 2)        # integer in {0, +-1, +-2}
        assert x.denominator == 1
        if (int(x) + 1) % 2 == 1 or t == u: r |= 1 << u
    rows.append(r)
piv = {}
for r in rows:
    x = r
    for p in sorted(piv, reverse=True):
        if x >> p & 1: x ^= piv[p]
    if x:
        p = x.bit_length() - 1
        for q in list(piv):
            if piv[q] >> p & 1: piv[q] ^= x
        piv[p] = x
assert len(piv) == 9
pv = sorted(piv)
def coords(r):
    c = 0
    for k, p in enumerate(pv):
        if r >> p & 1: r ^= piv[p]; c |= 1 << k
    assert r == 0
    return c
labT = [coords(r) for r in rows]
assert len(set(labT)) == 112
open(os.path.join(HERE, 'labT.txt'), 'w').write(' '.join(map(str, labT)) + '\n')
P = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../centre-floor-window') + '/')
L = open(P + 'pool_big.txt').read().split('\n')
n, Pn, F = map(int, L[0].split())
sets = [list(map(int, L[2 + i].split())) for i in range(Pn)]
mu = [float(x) for x in open(P + 'runBig/lpmwu_mu.txt').read().split()]
assert len(mu) == Pn
def det(rows):
    A = [[Fr(x) for x in r] for r in rows]; d = Fr(1)
    for c in range(8):
        p = next((r for r in range(c, 8) if A[r][c] != 0), None)
        if p is None: return Fr(0)
        if p != c: A[c], A[p] = A[p], A[c]; d = -d
        d *= A[c][c]
        for r in range(c + 1, 8):
            f = A[r][c] / A[c][c]
            if f: A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return d
ip = lambda a, b: sum(x * y for x, y in zip(a, b))
bad = 0; nfr = 0
out = []
for k, s in enumerate(sets):
    ms = [mp[t] for t in s]
    if len(set(ms)) != 8 or det([mine[t] for t in ms]) == 0: bad += 1
    if all(ip(mine[a], mine[b]) == 0 for i, a in enumerate(ms) for b in ms[i + 1:]): nfr += 1
    out.append(ms)
print('pool sets', Pn, 'non-bases', bad, 'orthogonal frames among them', nfr, '(first F =', F, ')')
c = [int(round(m * 2 ** 30)) for m in mu]
with open(os.path.join(HERE, 'pool_mine.txt'), 'w') as f:
    f.write('%d\n' % Pn)
    for ms, ci in zip(out, c): f.write('%d %s\n' % (ci, ' '.join(map(str, ms))))
print('sum c', sum(c))
