#!/usr/bin/env python3
"""v_thmA_raw.py (v-e3-slice6): the corrected U_3 program on FRESH random complex arrays, f = 1..6 digit columns.
Each register holds an array indexed by (x_1..x_f) in (F_2^3)^f; t_z acts as the diagonal i^(sum_c <z,x_c>) (the
X-basis form of the letter), X^a = t_a^2; u_z = ((1-i) + (1+i) X^z)/2 after each paid move (copy + shift + add);
outputs times prod_j ((1+i) + (1-i) X^ej)/2. Target: each output array = K a_w, K = t_e1 t_e2 t_e3. Registers: W arrays
plus M scratch (initially 0), moves are Lemma 1's conjugated rank-one updates."""
import numpy as np, time, sys
exec(open('v_thmA.py').read().split('def run(psi')[0])
rng = np.random.default_rng(int(time.time() * 1e6) % (2 ** 32))
for f in range(1, 7):
    Nf = 8 ** f
    xs = np.indices((8,) * f).reshape(f, -1)            # columns' digits
    def tdiag(z):
        e = sum(np.array([bin(z & x).count('1') % 2 for x in range(8)])[xs[c]] for c in range(f))
        return 1j ** (e % 4)
    a = rng.standard_normal((W, Nf)) + 1j * rng.standard_normal((W, Nf))
    R = np.zeros((W + M, Nf), complex); R[:W] = a
    for i in range(M):
        z = dirs[i]; t = tdiag(z); X = t * t; uz = ((1 - 1j) + (1 + 1j) * X) / 2
        y = V[:, i] @ R                                  # v_i^T R  (the register the move acts on, after mixing)
        R = R + np.outer(U[:, i], (t * uz - 1) * y)
    out = R[:W]
    for j in range(m):
        X = tdiag(1 << j) ** 2; out = out * ((1 + 1j) + (1 - 1j) * X) / 2
    K = tdiag(1) * tdiag(2) * tdiag(4)
    err = np.abs(out - K * a).max() / np.abs(a).max()
    print("f=%d: arrays of length %d, relative max error %.3g" % (f, Nf, err))
