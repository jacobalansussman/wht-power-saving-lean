"""symp.py extended: middle single-term contents whose term t sits at input distance d (d = 1: O1 output, O F O O F;
d = 2: O2 output, O O F O F; d = 3 for completeness).  M = sigma(t) L_Z at distance d, N' any complement of M."""
import sys, random
sys.path.insert(0, __import__('os').path.dirname(__file__))
from symp import *
from collections import Counter
rng = random.Random(5)
L1 = cell(G.K)
gens = [sigma(G.D[z]) for z in range(1, 16)] + [sigma(G.rand_engine(rng)) for _ in range(30)]
seen = {LZ}; fr = [LZ]
while fr:
    nf = []
    for L in fr:
        for S in gens:
            L2 = img(S, L)
            if L2 not in seen: seen.add(L2); nf.append(L2)
    fr = nf
Lag = list(seen)
for d in (1, 2, 3):
    Ms = [L for L in Lag if 4 - dimcap(L, LZ) == d]
    for name, els in (('Z7', z7()[0]), ('Z2^3', z2cube()[0])):
        C = [sigma(c) for c in els]
        cnt = Counter(); cex = 0; minfull = 99
        for M in Ms:
            CM = [img(c, M) for c in C]
            for Np in Lag:
                if len(M & Np) != 1: continue
                S = [g for g in range(len(C)) if img(C[g], Np) == Np]
                junk = [g for g in range(len(C)) if g not in S]
                full = [g for g in junk if len(CM[g] & Np) == 1]
                cnt[(len(S), len(full))] += 1
                if junk:
                    minfull = min(minfull, len(full))
                    if not full: cex += 1
        print(f'd={d} {name}: {len(Ms)} M x complements = {sum(cnt.values())} pairs; (|S|, #full-distance junk): {dict(sorted(cnt.items()))}; min full junk when S != H: {minfull}; R\' condition fails: {cex}', flush=True)
