#!/usr/bin/env python3
"""t_lift4.py (rm2-phase): an explicit frame family at block 4, built and checked here, as a test in C[Sp(8,2)].
A lift gives each of the 2295 frames U a 4-dim subspace Y_U of Q^8 spanned by integer rows reducing mod 2 to a basis
of U.  Opposite frames then have delta = 4 (mod-2 rank 8), every pair at distance d has delta >= d, and the bound is
abar_4/abar_1 = 4/(1 + E/34425), E = total edge excess.  Built greedily (BFS order, sampled sign patterns), then a few
improvement sweeps; final E recomputed with exact integer ranks; abar_4 = 4 re-checked on a random sample of opposite
pairs.  Usage: t_lift4.py [seed] [sweeps] [cands]."""
import sys, random, time, json
import numpy as np
from spg import Sp, lagrangians
from t_avg_lib import rank_int
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
sweeps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
cands = int(sys.argv[3]) if len(sys.argv) > 3 else 24
rng = random.Random(seed); t0 = time.time()
m = 4; sp = Sp(m); n = 2 * m
frames, idx = lagrangians(sp); F = len(frames)
mask = [sum(1 << u for u in U) for U in frames]
nb = [[] for _ in range(F)]
for i in range(F):
    mi = mask[i]
    for j in range(i + 1, F):
        if bin(mi & mask[j]).count("1") == 8: nb[i].append(j); nb[j].append(i)
edges = sum(len(x) for x in nb) // 2
print(f"block 4: {F} frames, {edges} edges, valency {set(len(x) for x in nb)}  ({time.time()-t0:.1f}s)"); sys.stdout.flush()
def basis_of(U):
    B = []; S = {0}
    vs = sorted(u for u in U if u); rng.shuffle(vs)
    for v in vs:
        if v not in S: B.append(v); S |= {a ^ v for a in S}
        if len(B) == m: break
    return B
def cand(U):
    B = basis_of(U)
    rows = []
    for v in B:
        rows.append([((v >> i) & 1) * rng.choice((1, -1)) for i in range(n)])
    return np.array(rows, dtype=np.int64)
def rk(A, Bm):
    return int(np.linalg.matrix_rank(np.vstack([A, Bm]).astype(float)))
Y = [None] * F
order = list(range(F))   # lagrangians() returns BFS order from L0
for k in order:
    best = None
    for _ in range(cands):
        c = cand(frames[k])
        e = sum(rk(c, Y[j]) - 5 for j in nb[k] if Y[j] is not None)
        if best is None or e < best[0]: best = (e, c)
        if e == 0: break
    Y[k] = best[1]
def total():
    return sum(rk(Y[i], Y[j]) - 5 for i in range(F) for j in nb[i] if j > i)
E = total(); print(f"greedy: E = {E}, ratio {4/(1+E/edges):.4f}  ({time.time()-t0:.1f}s)"); sys.stdout.flush()
for s in range(sweeps):
    perm = list(range(F)); rng.shuffle(perm)
    for k in perm:
        cur = sum(rk(Y[k], Y[j]) - 5 for j in nb[k])
        if cur == 0: continue
        for _ in range(cands):
            c = cand(frames[k]); e = sum(rk(c, Y[j]) - 5 for j in nb[k])
            if e < cur: Y[k] = c; cur = e
            if cur == 0: break
    E = total(); print(f"sweep {s+1}: E = {E}, ratio {4/(1+E/edges):.4f}  ({time.time()-t0:.1f}s)"); sys.stdout.flush()
# exact final check
Yl = [Y[k].tolist() for k in range(F)]
ok_red = True
for k in range(F):
    red = [sum((abs(x) % 2) << i for i, x in enumerate(r)) for r in Yl[k]]
    S = {0}
    for v in red: S |= {a ^ v for a in S}
    ok_red &= (frozenset(S) == frames[k])
Eex = sum(rank_int(Yl[i] + Yl[j]) - 5 for i in range(F) for j in nb[i] if j > i)
opp = [(i, j) for i in range(F) for j in range(F) if bin(mask[i] & mask[j]).count("1") == 1]
smp = rng.sample(opp, 3000)
ok_opp = all(rank_int(Yl[i] + Yl[j]) == 8 for i, j in smp)
from fractions import Fraction as Fr
ab1 = 1 + Fr(Eex, edges); ratio = Fr(4) / ab1
print(f"EXACT: rows reduce to their frames {ok_red}; total edge excess E = {Eex}; abar_1 = {ab1}; opposite pairs delta = 4 "
      f"(sample of 3000 of {len(opp)}; forced by the mod-2 rank anyway) {ok_opp}; ratio abar_4/abar_1 = {ratio} = {float(ratio):.6f}")
print(f"   so in C[Sp(8,2)] (Lemma N's block) every unit has M >= {float(ratio):.4f} W; 4 in 15 needs > 3.75: refused {ratio * 4 > 15}")
json.dump({"m": 4, "E": Eex, "frames_pauli_labels": [sorted(U) for U in frames], "rows": Yl}, open(f"lift_m4_s{seed}.json", "w"))
print(f"done ({time.time()-t0:.1f}s)")
