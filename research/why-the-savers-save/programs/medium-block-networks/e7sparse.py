#!/usr/bin/env python3
"""cx-medium / e7sparse.py (own code; numpy from our pylib; python3 -I).  Is there an E7 decoder below the one-pair
break-even (42 per call) in the hub's accounting?  E7's 63 lines = the nonzero vectors e of a 6-dim symplectic F_2^6
(label u = 1 + e in F_2^7; u_a . u_b = 1 + w(e_a, e_b)); a pair must be exact iff it commutes (w = 0).
A decoder is a gather matrix G (t x 63) and scatter rows; in the hub's accounting a total costs dim span{u_b : b in
supp G_k} (its copy descends to frame 0), whatever it scatters to.  A gather matrix G admits a valid scatter iff for
every port a, column G[:, a] is NOT in the span of the columns of the 30 ports commuting with a (then a scatter row
r_a with r_a . G[:, b] = delta_ab on a's commuting set exists, by linear algebra).
cliquefloor.py: every decoder costs >= 35 per call; its bound is met only by totals gathering the 15 ports of 1 + V,
V a 4-dim subspace (cost 5).  Here: random searches over supports of t totals inside 4-dim (cost 5) and 3-dim (cost 4)
subspaces, generic real values (and +-1 values), testing the column condition exactly by rank (numpy, tolerance), for
total cost 35 .. 41."""
import sys, os, json, itertools, random
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from designs_cx import design, labels, dot

L = design('E7'); labs, k, B = labels(L); n = len(L)
ONE = (1 << k) - 1
E = [u ^ ONE for u in labs]                          # even parts e = u + 1 (6-dim symplectic space inside F_2^7)
comm = [[b for b in range(n) if b != a and B[a][b] == 1] for a in range(n)]
assert all(len(c) == 30 for c in comm)

def span(vs):
    S = {0}
    for v in vs: S |= {x ^ v for x in S}
    return S

evens = sorted(set(E))
def rand_subspace(d, rng):
    while True:
        gens = rng.sample(evens, d)
        S = span(gens)
        if len(S) == 1 << d: return S
def ports_in(S): return [a for a in range(n) if E[a] in S]

def valid(G, tol=1e-9):
    bad = 0
    for a in range(n):
        Ca = comm[a]
        r0 = np.linalg.matrix_rank(G[:, Ca], tol=tol) if Ca else 0
        r1 = np.linalg.matrix_rank(G[:, Ca + [a]], tol=tol)
        if r1 == r0: bad += 1
    return bad

def cost_of(supp):
    labsS = [labs[a] for a in supp]
    b = []
    for v in labsS:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v)
    return len(b)

rng = random.Random(15); nrng = np.random.default_rng(15)
results = []
best = None
plans = [(7, 0), (8, 0), (6, 2), (5, 4), (7, 1), (4, 5), (3, 6), (2, 7)]    # (#cost-5 totals, #cost-4 totals)
TR = int(os.environ.get('TRIALS', '400'))
for n5, n4 in plans:
    total_cost = 5 * n5 + 4 * n4
    if total_cost > 41: continue
    fewest = None
    for trial in range(TR):
        # greedy covers: each new total takes, among 300 random subspaces of its size, the one covering most new ports
        sup, cover = [], set()
        for d in [4] * n5 + [3] * n4:
            cands = [ports_in(rand_subspace(d, rng)) for _ in range(300)]
            pick = max(cands, key=lambda s: (len(set(s) - cover), rng.random()))
            sup.append(pick); cover |= set(pick)
        cover = set().union(*map(set, sup))
        if len(cover) < n: continue
        G = np.zeros((len(sup), n))
        for i, s in enumerate(sup):
            G[i, s] = nrng.normal(size=len(s))
        b = valid(G)
        if fewest is None or b < fewest: fewest = b
        if b == 0:
            best = dict(n5=n5, n4=n4, cost=total_cost, supports=sup); break
    results.append(dict(n5=n5, n4=n4, cost=total_cost, fewest_bad_ports=fewest))
    print(f"{n5} totals of cost 5 + {n4} of cost 4 (cost {total_cost}): fewest ports failing the column condition over"
          f" {TR} random covers = {fewest}", flush=True)
    if best: break
print('valid decoder below 42 found:', bool(best), best and best['cost'])
json.dump(dict(results=results, found=best), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'e7sparse.json'), 'w'), indent=1, default=str)
