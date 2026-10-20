#!/usr/bin/env python3
"""ctr-window / worstW.py BASIS...: describe the gather set L(W) of a label subspace W (basis given as ints): number of
lines, rank of their roots over Q, and how many pairs are orthogonal / at 60 degrees; whether the roots form a closed
root subsystem (reflections) within E8, and whether it avoids the removed frame."""
import sys, os
from fractions import Fraction
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from designs_cx import design, labels, ip, roots_E8, lines
L = design('P112'); labs, k, _ = labels(L); lab2 = {u: p for p, u in enumerate(labs)}
B = [int(x) for x in sys.argv[1:]]
els = {0}
for v in B: els |= {x ^ v for x in els}
S = sorted(lab2[x] for x in els if x in lab2)
def rankQ(vecs):
    rows = [[Fraction(x) for x in v] for v in vecs]; r = 0
    for c in range(8):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]; rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r
roots = {tuple(v) for t in S for v in (L[t], tuple(-x for x in L[t]))}
closed = all(tuple(a[i] - 2 * ip(a, b) // ip(b, b) * b[i] for i in range(8)) in roots for a in roots for b in roots)
E8lines = lines(roots_E8()); live = set(map(tuple, L))
span_r = rankQ([L[t] for t in S])
# lines of E8 (incl. removed frame) in the Q-span closure that are not in S
orth = sum(1 for i in S for j in S if i < j and ip(L[i], L[j]) == 0)
print('dim W', len(B), '| |L(W)|', len(S), '| root rank', span_r, '| closed root subsystem', closed, '| orthogonal pairs', orth, 'of', len(S) * (len(S) - 1) // 2)
