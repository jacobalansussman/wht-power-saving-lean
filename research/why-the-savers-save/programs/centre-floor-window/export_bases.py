#!/usr/bin/env python3
"""ctr-window / export_bases.py N SEED: write N random 8-subsets of the P112 lines whose roots form a basis of R^8
(exact rank over Q), in the same format as p112.txt, to bases_N_SEED.txt. In a rank-8 (Gram up to gauge) realisation
each such S has N[S,S] = D1 G[S,S] D2 of rank 8, so every single-layer decoder needs >= 8 terms meeting S on both sides."""
import sys, os, random
from fractions import Fraction
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from designs_cx import design, labels, ip
L = design('P112'); n = len(L); labs, k, _ = labels(L)
N = int(sys.argv[1]); seed = int(sys.argv[2])
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
rng = random.Random(seed); out = []
while len(out) < N:
    S = sorted(rng.sample(range(n), 8))
    if rankQ([L[s] for s in S]) == 8: out.append(S)
with open('bases_%d_%d.txt' % (N, seed), 'w') as f:
    f.write('%d %d\n' % (n, N)); f.write(' '.join(map(str, labs)) + '\n')
    for S in out: f.write(' '.join(map(str, S)) + '\n')
print('wrote', len(out))
