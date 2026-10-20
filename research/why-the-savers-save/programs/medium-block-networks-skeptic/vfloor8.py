#!/usr/bin/env python3
"""v-cx-medium / vfloor8.py (own code, stdlib; python3 -I).  Exact clique floor for E8's 56-line hyperplane section
(k = 8): enumerate EVERY subspace of the 8-dim label space (417,199, by reduced row echelon form), count the maximum
cliques met by the ports whose label coordinates lie in it, and take min dim/T.  Label coordinates = coordinates of the
rows of B in a basis of its row space (ranks preserved).  Calibrated first on E7 (k = 7): must give 35."""
import sys, itertools
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from vclique import design, gram_B, rowint, add_basis, max_cliques
from vsub2 import section
from fractions import Fraction as Fr

def coords_of(L):
    n = len(L); B = gram_B(L); rows = [rowint(B[i]) for i in range(n)]
    basis = {}
    for r in rows: add_basis(basis, r)
    bv = list(basis.values())
    piv = {b.bit_length() - 1: (i, b) for i, b in enumerate(bv)}
    def co(v):
        c = 0
        while v:
            p = v.bit_length() - 1; i, b = piv[p]; v ^= b; c ^= 1 << i
        return c
    return [co(r) for r in rows], len(bv), B

def rref_subspaces(k):
    """yield (d, list of basis ints) for every subspace of F_2^k, via RREF with pivot columns chosen"""
    for d in range(1, k + 1):
        for piv in itertools.combinations(range(k), d):
            # free positions: for row i (pivot piv[i]), columns > piv[i] that are not pivots
            free = [[c for c in range(piv[i] + 1, k) if c not in piv] for i in range(d)]
            nf = sum(len(f) for f in free)
            for bits in range(1 << nf):
                rows = []; t = bits
                for i in range(d):
                    v = 1 << piv[i]
                    for c in free[i]:
                        if t & 1: v |= 1 << c
                        t >>= 1
                    rows.append(v)
                yield d, rows

def floor(name, L):
    C, k, B = coords_of(L); n = len(L)
    adj = [{j for j in range(n) if j != i and B[i][j] == 1} for i in range(n)]
    om, cl = max_cliques(n, adj)
    of = [0] * n
    for ci, c in enumerate(cl):
        for p in c: of[p] |= 1 << ci
    lab = {}
    for i, c in enumerate(C): lab.setdefault(c, []).append(i)
    best = None; count = 0
    for d, rows in rref_subspaces(k):
        count += 1
        S = [0]
        for v in rows: S = S + [s ^ v for s in S]
        m = 0
        for s in S:
            for p in lab.get(s, ()): m |= of[p]
        T = bin(m).count('1')
        if T and (best is None or Fr(d, T) < best[0]): best = (Fr(d, T), d, T)
    f = om * len(cl) * best[0]
    print(f"{name}: n={n} k={k} omega={om} cliques={len(cl)} subspaces={count} best dim/T={best[1]}/{best[2]} floor={float(f):.3f} per call"
          f" ({float(f)/n:.3f} per port); one-pair break-even {2*n/3:.2f}", flush=True)
    return f

if __name__ == '__main__':
    floor('E7 (calibration, must be 35)', design('E7'))
    floor('E8-56', section(design('E8'), 56, 4))
