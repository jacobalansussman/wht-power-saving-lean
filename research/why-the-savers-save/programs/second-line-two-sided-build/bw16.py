#!/usr/bin/env python3
"""an-build: an-synth proposal 2 screen on the Barnes-Wall lattice BW16 (own code, stdlib).  Minimal vectors at norm 8:
(+-2,+-2,0^14) on any two coordinates (480) and (+-1^8, 0^8) on the 30 weight-8 words of RM(1,4) with an even number of
minus signs (3840); 4320 vectors, 2160 lines.  Screen (an-synth 4.2): centres t = rank A, tightness (A^2 = (n/t) A <=>
sum f f^T = (n/t) I), junk per port (non-orthogonal other lines), a junk clique (lower bound on the window k), copy ratio
bound t*k/n, against the pass marks: copy ratio <= 0.35, window <= 12, junk per port <= 120."""
import itertools, json
from fractions import Fraction as F
from collections import Counter
# RM(1,4): affine functions on F_2^4
pts = list(itertools.product([0, 1], repeat=4))
words = set()
for a in itertools.product([0, 1], repeat=4):
    for c in (0, 1):
        w = tuple((sum(ai * pi for ai, pi in zip(a, p)) + c) % 2 for p in pts)
        words.add(w)
octads = [w for w in words if sum(w) == 8]
assert len(octads) == 30
vecs = []
for i, j in itertools.combinations(range(16), 2):
    for si, sj in itertools.product([2, -2], repeat=2):
        v = [0] * 16; v[i] = si; v[j] = sj; vecs.append(tuple(v))
for w in octads:
    sup = [k for k in range(16) if w[k]]
    for signs in itertools.product([1, -1], repeat=8):
        if signs.count(-1) % 2 == 0:
            v = [0] * 16
            for k, s in zip(sup, signs): v[k] = s
            vecs.append(tuple(v))
assert len(vecs) == 4320 and all(sum(x * x for x in v) == 8 for v in vecs)
# lines: one representative per +-pair
seen, lines = set(), []
for v in vecs:
    m = tuple(-x for x in v)
    if m not in seen:
        seen.add(v); lines.append(v)
n = len(lines); assert n == 2160
ip = lambda a, b: sum(x * y for x, y in zip(a, b))
# closure check: is the set of minimal vectors closed under negation and a lattice's minimal set? (count inner products)
dist = Counter(ip(lines[0], w) for w in lines[1:])
junk = sum(c for k, c in dist.items() if k != 0)
# tightness: sum over lines of f f^T is a multiple of I  (frame operator), with f = v / sqrt(8)
S = [[sum(v[i] * v[j] for v in lines) for j in range(16)] for i in range(16)]
tight = all(S[i][j] == (S[0][0] if i == j else 0) for i in range(16) for j in range(16))
t = 16                                              # the lines span R^16 (rank of A = rank of the frame operator)
# a junk clique: lines 2e_0 + 2e_j, j = 1..15, pairwise inner product 4 (normalised 1/2)
clique = [tuple([2] + [2 if k == j else 0 for k in range(1, 16)]) for j in range(1, 16)]
assert all(ip(a, b) != 0 for a, b in itertools.combinations(clique, 2))
# greedy extension of the clique inside the junk graph
cl = list(clique)
for w in lines:
    if w not in cl and all(ip(w, c) != 0 for c in cl): cl.append(w)
omega_lb = len(cl)
out = dict(n_lines=n, centres_t=t, n_over_t=n / t, tight=tight, inner_products_with_one_line=dict(sorted(dist.items())),
           junk_per_port=junk, junk_clique_lower_bound=omega_lb, window_lower_bound=omega_lb,
           copy_ratio_lower_bound=t * omega_lb / n,
           passes=dict(copy_ratio=t * omega_lb / n <= 0.35, window=omega_lb <= 12, junk=junk <= 120))
print(json.dumps(out, indent=1))
json.dump(out, open('bw16.json', 'w'), indent=1)
