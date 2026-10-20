#!/usr/bin/env python3
"""s14_fixedmove_check.py (t2-block3): calibration of PROOF.md S4 for fixed moves.
c = 3-cycle (commutes with K, so c' = c^-1); a fixed array y restricted to N_k = { y in N : c' y = zeta^k y }; the raw
input register is tau-fixed and dir7 commutes with c.  The potential gain of the fixed move on these inputs is
dim((N + rho(dir7) N_k)/N); S4 predicts it is at most sigma_{-k}, sigma = split of rho(c) on (N + rho(dir7) N)/N."""
import numpy as np
from lib3 import *
def coordperm(p):
    def f(x):
        y = 0
        for j in range(3):
            if (x >> j) & 1: y |= 1 << p[j]
        return y
    return perm_mat(f)
c = coordperm([1, 2, 0]); cp = Kinv @ np.linalg.inv(c) @ K
assert np.allclose(cp, np.linalg.inv(c))
Nb = np.zeros((560, 280), dtype=complex); Nb[NIDX, np.arange(280)] = 1
Rcp = rho(cp); Rm = rho(dirz(7)); Rc = rho(c)
z = np.exp(2j * np.pi / 3)
NM = colspace(np.hstack([Nb, Rm @ Nb]))
sig = []
for k in range(3):
    P = sum(z ** (-k * s) * np.linalg.matrix_power(Rc, s) for s in range(3)) / 3
    sig.append(colspace(P @ NM).shape[1] - colspace(P @ Nb).shape[1])
print("sigma (c-eigenvalue zeta^k on (N+mN)/N), k = 0,1,2:", sig)
for k in range(3):
    P = sum(z ** (-k * s) * np.linalg.matrix_power(Rcp, s) for s in range(3)) / 3
    Nk = colspace(P @ Nb)
    gain = colspace(np.hstack([Nb, Rm @ Nk])).shape[1] - 280
    print("k = %d: dim N_k = %d, gain of the fixed move on N_k-inputs = %d, predicted bound sigma_{-k} = %d" % (k, Nk.shape[1], gain, sig[(-k) % 3]))
