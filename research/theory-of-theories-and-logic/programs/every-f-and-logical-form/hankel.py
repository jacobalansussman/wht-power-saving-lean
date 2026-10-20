"""hankel.py -- ml-logic (meta3). Which letters are finite-state across digit columns?
A family of diagonal letters d_f(w) (w = the f column digits) is a weighted automaton (MPO of bounded bond dimension)
iff the Hankel matrix H[u, v] = d(uv) has finite rank (Fliess / Carlyle-Paz).  We compare
  (1) a column-uniform phase  d(w) = prod_t i^{w_t}                    -> rank 1
  (2) a neighbour phase (CZ chain) d(w) = prod_t (-1)^{w_t w_{t+1}}    -> rank 2
  (3) a Cooley-Tukey twiddle  d(w) = exp(2 pi i j * val(w)), val(w) = sum_t w_t 2^-t (digits read most significant
      first, block m = 1), j = 1                                       -> rank grows with the length: not finite-state.
Hankel blocks over all u, v of length <= L (2^(L+1)-1 rows/cols), numerical rank with tolerance 1e-9.
Usage: python3 -I hankel.py   (tiny)
"""
import sys
import numpy as np, itertools, cmath

def words(L):
    out = [()]
    for l in range(1, L + 1): out += list(itertools.product((0, 1), repeat=l))
    return out

def d1(w): return 1j ** sum(w)
def d2(w): return (-1) ** sum(w[t] * w[t + 1] for t in range(len(w) - 1))
def d3(w): return cmath.exp(2j * cmath.pi * sum(b * 2.0 ** -(t + 1) for t, b in enumerate(w)))

for L in (2, 3, 4, 5, 6):
    W = words(L); ranks = []
    for d in (d1, d2, d3):
        H = np.array([[d(u + v) for v in W] for u in W])
        s = np.linalg.svd(H, compute_uv=False)
        ranks.append(int((s > 1e-9 * s[0]).sum()))
    print(f"L = {L}: Hankel ranks  uniform phase {ranks[0]},  CZ chain {ranks[1]},  twiddle {ranks[2]}", flush=True)
