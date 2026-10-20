#!/usr/bin/env python3
"""v-cx-medium / vsub3.py (own code; numpy; python3 -I).  E8's 56-line hyperplane section (k = 8, omega = 4; one pair ->
block 22): the fast facts.  (1) The clique floor can never exceed omega * k (take W = the whole label space: dim k,
meets every clique), so here it is <= 32 < 37.33 = 2n/3: cx-medium's floor tool cannot exclude this template point.
(2) Coordinate / symmetric frame decoders: cost of the coordinate decoder over many orthonormal bases; and the symmetric
bound: each total gathers >= n - (max lines in a hyperplane of R^8), whose labels need dim >= ..."""
import sys, json, itertools
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from vclique import design, gram_B, rowint, add_basis, max_cliques
from vsub2 import section
L = section(design('E8'), 56, 4); n = len(L)
B = gram_B(L); rows = [rowint(B[i]) for i in range(n)]
def rk(S):
    b = {}
    for p in S: add_basis(b, rows[p])
    return len(b)
k = rk(range(n))
adj = [{j for j in range(n) if j != i and B[i][j] == 1} for i in range(n)]
om, cl = max_cliques(n, adj)
junk = sorted(set(sum(1 for j in range(n) if B[i][j] == 0) for i in range(n)))
orth = sorted(set(len(adj[i]) for i in range(n)))
sig = sorted(set(rk([j for j in range(n) if B[i][j] == 0]) for i in range(n)))
R = np.array([[float(x) for x in l] for l in L])
G = R @ R.T
print(f"n={n} k={k} omega={om} maxcliques={len(cl)} junk/port {junk} orthogonal partners/port {orth} sigma {sig}; "
      f"span dim {np.linalg.matrix_rank(R)}; floor <= omega*k = {om*k}; break-even (one pair) {2*n/3:.2f}; block {1+3*(k-1)}")
# max lines of the section in a hyperplane of R^8: hyperplanes spanned by 7 independent lines (exhaustive over the
# hyperplanes containing at least 7 independent lines, via all 7-subsets is too many; use: the lines orthogonal to f
# for f = each line, each sum/difference of two lines, random f)
best_perp = 0
cands = [R[i] for i in range(n)] + [R[i] + s * R[j] for i in range(n) for j in range(i + 1, n) for s in (1, -1)]
rng = np.random.default_rng(1)
cands += list(rng.normal(size=(2000, 8)))
for f in cands:
    c = int(sum(1 for i in range(n) if abs(R[i] @ f) < 1e-9))
    best_perp = max(best_perp, c)
print('max lines orthogonal to one f (sampled/structured):', best_perp, '-> each symmetric total gathers >=', n - best_perp)
# coordinate decoders: random orthonormal bases + structured ones (the E8 coordinate basis and Hadamard basis)
Hd = np.array([[1.0]])
for _ in range(3): Hd = np.block([[Hd, Hd], [Hd, -Hd]])
bases = [np.eye(8), Hd / np.sqrt(8)]
for t in range(200):
    A, _ = np.linalg.qr(rng.normal(size=(8, 8))); bases.append(A.T)
best = None
for Fm in bases:
    Gm = Fm @ R.T
    M = Gm.T @ Gm / 2
    assert all(abs(M[a, b] - (1.0 if a == b else 0.0)) < 1e-9 for a in range(n) for b in range(n) if B[a][b] == 1)
    cost = [rk([b for b in range(n) if abs(Gm[q, b]) > 1e-9]) for q in range(8)]
    if best is None or sum(cost) < best[0]: best = (sum(cost), cost)
print('best coordinate decoder over', len(bases), 'orthonormal bases:', best)
# label rank needed for g gathered labels: max number of section labels in a d-dim label-spanned subspace (greedy upper
# estimate by closures from random seeds)
