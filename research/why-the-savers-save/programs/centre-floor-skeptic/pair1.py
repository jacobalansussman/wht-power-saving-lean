#!/usr/bin/env python3
"""v-centre / pair1.py: own check that pairing-rank-1 terms can meet ~600 of the 1260 frames on both sides
(so frame counting gives no useful floor in the pairing-rank accounting).  R = an 8-line set, A = the lines whose
B-column is constant on R; then rank_F2 B[R,A] <= 1.  Count frames meeting both R and A."""
import json, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'data.json')))
L = [tuple(v) for v in D['lines']]; frames = [set(f) for f in D['frames']]
ip = lambda a, b: sum(x * y for x, y in zip(a, b))
B = [[1 if (t == s or ip(L[t], L[s]) == 0) else 0 for s in range(112)] for t in range(112)]
def rank2(rows):
    b = {}
    for x in rows:
        for p in sorted(b, reverse=True):
            if x >> p & 1: x ^= b[p]
        if x: b[x.bit_length() - 1] = x
    return len(b)
best = (0, None)
def test(R):
    A = [a for a in range(112) if len({B[r][a] for r in R}) == 1]
    rk = rank2([sum(B[r][a] << i for i, a in enumerate(A)) for r in R])
    Rs, As = set(R), set(A)
    hit = sum(1 for f in frames if f & Rs and f & As)
    return hit, rk, len(A)
for f in D['frames']:
    h, rk, na = test(f)
    if h > best[0]: best = (h, ('frame', f, rk, na))
random.seed(1)
for _ in range(3000):
    R = random.sample(range(112), random.choice([2, 3, 4, 6, 8]))
    h, rk, na = test(R)
    if rk <= 1 and h > best[0]: best = (h, ('rand', R, rk, na))
print('best frames met on both sides by a rank<=1 term:', best[0], best[1], ' -> counting bound <=', round(8 * 1260 / best[0], 2))
