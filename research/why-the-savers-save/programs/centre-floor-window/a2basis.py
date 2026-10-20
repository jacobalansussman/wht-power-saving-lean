#!/usr/bin/env python3
"""ctr-window / a2basis.py (own code; python3 -I). For the best terms found by a2clique.py (pairs of label subspaces
W1 >= W2^perp), report |R| = #lines with label in W1, |A|, the pairing rank, the exact-clique hit count, and the
fraction of random spanning 8-sets (bases of R^8 drawn uniformly from 8-subsets of the 112 lines, rejection-sampled)
that meet both R and A. In a rank-8 (Gram, up to gauge) realisation every basis S needs >= 8 terms meeting S on both
sides, so a term ratio fraction/rank above 1/8 means basis counting cannot reach 64 either."""
import sys, os, json, random
from fractions import Fraction
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from designs_cx import design, labels, ip
L = design('P112'); n = len(L); labs, k, _ = labels(L)
lab2port = {u: p for p, u in enumerate(labs)}
def dot(x, y): return bin(x & y).count('1') & 1
def span(vs):
    els = {0}
    for v in vs:
        if v not in els: els |= {x ^ v for x in els}
    return els
def rankQ(vecs):
    rows = [[Fraction(x) for x in v] for v in vecs]; r = 0; cols = len(rows[0]) if rows else 0
    for c in range(cols):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]; rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), sys.argv[1] if len(sys.argv) > 1 else 'a2clique.json')))
rng = random.Random(7)
bases = []
while len(bases) < 3000:
    S = rng.sample(range(n), 8)
    if rankQ([L[s] for s in S]) == 8: bases.append(set(S))
for r, v in sorted(d['A2_search'].items(), key=lambda kv: int(kv[0])):
    W1 = span([int(x, 16) for x in v['W1']]); W2 = span([int(x, 16) for x in v['W2']])
    R = {lab2port[x] for x in W1 if x in lab2port}; A = {lab2port[x] for x in W2 if x in lab2port}
    frac = sum(1 for S in bases if S & R and S & A) / len(bases)
    print('rank %s: |R| = %d, |A| = %d, clique hits %d, basis-hit fraction %.3f, ratio/(1/8) = %.2f'
          % (r, len(R), len(A), v['hits'], frac, frac / int(r) * 8))
