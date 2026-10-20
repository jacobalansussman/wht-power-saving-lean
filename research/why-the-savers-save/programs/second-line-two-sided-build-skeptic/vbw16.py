#!/usr/bin/env python3
"""v-an-build: BW16 screen, own code, independent construction.  Enumerate ALL norm-8 vectors of Z^16 of the three
possible shapes (2,2,0^14), (2,1^4,0^11), (1^8,0^8) with signs, keep x with (x mod 2) in RM(1,4) and sum(x) = 0 mod 4
(the standard definition of BW16), and screen: lines, inner-product distribution, junk per port, rank of the Gram
(centres t), tightness, a large mutually-junk clique (greedy + randomised restarts).  Then the Leech clique {4e1+4ej}."""
import itertools, random, json
from fractions import Fraction as F
from collections import Counter
pts = list(itertools.product([0, 1], repeat=4))
RM = set()
for a in itertools.product([0, 1], repeat=4):
    for c in (0, 1):
        RM.add(tuple((sum(x * y for x, y in zip(a, p)) + c) % 2 for p in pts))
def ok(x):
    return tuple(t % 2 for t in x) in RM and sum(x) % 4 == 0
V = set()
for sup in itertools.combinations(range(16), 2):
    for s in itertools.product([2, -2], repeat=2):
        x = [0] * 16
        for i, t in zip(sup, s): x[i] = t
        if ok(x): V.add(tuple(x))
for i in range(16):
    for sup in itertools.combinations([j for j in range(16) if j != i], 4):
        for s0 in (2, -2):
            for s in itertools.product([1, -1], repeat=4):
                x = [0] * 16; x[i] = s0
                for j, t in zip(sup, s): x[j] = t
                if ok(x): V.add(tuple(x))
n24 = 0
for sup in itertools.combinations(range(16), 8):
    if tuple(1 if k in sup else 0 for k in range(16)) not in RM: continue
    for s in itertools.product([1, -1], repeat=8):
        x = [0] * 16
        for j, t in zip(sup, s): x[j] = t
        if ok(x): V.add(tuple(x))
V = sorted(V)
lines = sorted(set(max(v, tuple(-t for t in v)) for v in V))
ip = lambda a, b: sum(x * y for x, y in zip(a, b))
print('norm-8 vectors in BW16:', len(V), ' lines:', len(lines), ' shapes:', Counter(tuple(sorted(Counter(abs(t) for t in v if t).items())) for v in V))
dist = Counter()
for L in lines[:50]:
    c = Counter(abs(ip(L, M)) for M in lines if M != L)
    dist[tuple(sorted(c.items()))] += 1
print('|ip| distribution from each of 50 lines (|ip| 0, 2, 4 of norm 8 = 0, 1/4, 1/2):', dict(dist))
# tightness: sum over lines of f f^T = (n/16) I  (frame), checked exactly
G = [[0] * 16 for _ in range(16)]
for L in lines:
    for i in range(16):
        if L[i]:
            for j in range(16): G[i][j] += L[i] * L[j]
diag = set(G[i][i] for i in range(16)); off = set(G[i][j] for i in range(16) for j in range(16) if i != j)
print('sum f f^T: diagonal', diag, 'off-diagonal', off, '=> tight, n/t =', F(len(lines), 16))
# clique of mutually non-orthogonal lines
adj = {i: set(j for j in range(len(lines)) if j != i and ip(lines[i], lines[j]) != 0) for i in range(len(lines))}
random.seed(1)
best = []
for trial in range(300):
    order = list(range(len(lines))); random.shuffle(order)
    cl = []
    cand = set(range(len(lines)))
    for i in order:
        if i in cand:
            cl.append(i); cand &= adj[i]
    if len(cl) > len(best): best = cl
assert all(ip(lines[a], lines[b]) != 0 for a, b in itertools.combinations(best, 2))
print('mutually-junk clique found:', len(best), '(window >= this; pass mark <= 12)')
# an explicit structured clique: all (1^8) lines on one octad with even minus count, first sign +: check
# Leech: 4e1 + 4ej, j = 2..24
leech = [tuple(4 if k in (0, j) else 0 for k in range(24)) for j in range(1, 24)]
print('Leech {4e1+4ej}: %d lines, pairwise ip %s' % (len(leech), set(ip(a, b) for a, b in itertools.combinations(leech, 2))))
json.dump(dict(vectors=len(V), lines=len(lines), ipdist={str(k): v for k, v in dist.items()}, clique=len(best),
               clique_lines=[lines[i] for i in best]), open('vbw16.json', 'w'))
