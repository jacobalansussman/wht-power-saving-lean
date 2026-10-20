#!/usr/bin/env python3
"""vpolyalt.py: deletion bound max|J|/|K| in the polynomial-basis placement for EVERY primitive polynomial, m=3..5
(robustness of S5(b), which is claimed for the polynomials of e3-rmbuild/w8)."""
import sys, itertools
sys.path.insert(0, '.')
import vfield
from vfield import setup, channel_data
for m in (3, 4, 5):
    res = []
    for p in range(1 << m, 1 << (m + 1)):
        if not (p & 1): continue
        vfield.PRIM[m] = p
        try:
            F, basis, C, sig, Sk, pairs = setup(m, 'poly')
        except AssertionError:
            continue
        N = len(pairs); best = 0
        for q in itertools.product([0, 1], repeat=N):
            if not any(q): continue
            K, J = channel_data(F, C, Sk, pairs, q)
            if sum(J): best = max(best, sum(J) / sum(K))
        res.append((bin(p), round(best, 6)))
    print(f"m={m}: poly-placement deletion bound per primitive polynomial: {res}", flush=True)
