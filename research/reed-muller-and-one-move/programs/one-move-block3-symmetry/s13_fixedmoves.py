#!/usr/bin/env python3
"""s13_fixedmoves.py (t2-block3): robustness of the fixed-move cap splits used by s5/s11/s12.
For T = <transposition (01)>, <3-cycle>, S3 (coordinate permutations) and <n> (the F21 element), collect every move
letter m = b dir z b' (b from a menu of products of <= 2 free generators, b' in {1, b^-1}) that commutes with T up to a
scalar, and the distinct traces of T's elements on (N + rho(m) N)/N."""
import itertools
import numpy as np
from lib3 import *

def coordperm(p):
    def f(x):
        y = 0
        for j in range(3):
            if (x >> j) & 1: y |= 1 << p[j]
        return y
    return perm_mat(f)
t01 = coordperm([1, 0, 2]); c3 = coordperm([1, 2, 0])
gens = [X(1), X(2), X(4), Z(1), Z(2), Z(4), phase_letter(1), phase_letter(2), phase_letter(4), phase_letter(7),
        phase_letter(3), phase_letter(0, 1)]
for i in range(3):
    for j in range(3):
        if i != j:
            gens.append(perm_mat(lambda x, i=i, j=j: x ^ (((x >> i) & 1) << j)))
menu = [I8] + gens + [g @ h for g in gens for h in gens]
Nb = np.zeros((560, 280), dtype=complex); Nb[NIDX, np.arange(280)] = 1
groups = {"<(01)>": [t01], "<3-cycle>": [c3], "S3": [t01, c3]}
cache = {}
def commutes(g, m):
    r = g @ m @ np.linalg.inv(g)
    k0 = np.argmax(np.abs(m.flatten()))
    sc = r.flatten()[k0] / m.flatten()[k0]
    return np.allclose(r, sc * m)
results = {k: {} for k in groups}
seen = set()
for z in range(1, 8):
    D = dirz(z)
    for b in menu:
        for m in (b @ D @ np.linalg.inv(b), b @ D):
            kk = key(m)
            if kk in seen: continue
            seen.add(kk)
            for gname, gl in groups.items():
                if all(commutes(g, m) for g in gl):
                    if kk not in cache:
                        NM = colspace(np.hstack([Nb, rho(m) @ Nb]))
                        cache[kk] = NM
                    NM = cache[kk]
                    pN = np.linalg.pinv(NM)
                    trs = tuple(int(round((np.trace(pN @ rho(g) @ NM) - np.trace(Nb.conj().T @ rho(g) @ Nb)).real)) for g in [t01, c3])
                    results[gname].setdefault(trs, 0)
                    results[gname][trs] += 1
for gname, d in results.items():
    print("%-10s fixed move letters found: %d; distinct (trace of (01), trace of 3-cycle) on (N+mN)/N: %s" % (
        gname, sum(d.values()), sorted(d.items())))
