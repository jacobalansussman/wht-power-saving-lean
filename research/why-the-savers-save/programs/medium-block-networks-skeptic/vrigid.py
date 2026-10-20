#!/usr/bin/env python3
"""v-cx-medium / vrigid.py (own code; numpy from our pylib; python3 -I).  Independent check of cx-medium's e7rigid claim:
at the root (coordinate) decoder of E7, the tangent space of {(Phi, Psi) in R^{7x63} x R^{7x63} : (Psi^T Phi)_ab = delta_ab
on every constrained entry (diagonal and orthogonal root pairs)} has dimension 111 = 49 (GL(7)) + 62 (column rescaling
Phi -> Phi D, Psi -> Psi D^-1, minus the scalar shared with GL(7)).  My E7 is the E8 roots orthogonal to (1/2)^8, with
the Hadamard basis of its span (cx-medium: E8 roots orthogonal to e1+e2, coordinate basis)."""
import sys, itertools
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from vclique import design, ip
L = design('E7'); n = len(L)
R = np.array([[float(x) for x in l] for l in L])
Hd = np.array([[1.0]])
for _ in range(3): Hd = np.block([[Hd, Hd], [Hd, -Hd]])
F = Hd[1:] / np.sqrt(8)                     # 7 orthonormal vectors orthogonal to (1,...,1)
Phi = F @ R.T                              # 7 x 63
Psi = Phi / 2.0
M = Psi.T @ Phi
G = R @ R.T
cons = [(a, b) for a in range(n) for b in range(n) if a == b or abs(G[a, b]) < 1e-9]
assert all(abs(M[a, b] - (1.0 if a == b else 0.0)) < 1e-12 for a, b in cons)
print('E7 lines', n, 'constrained entries', len(cons), 'junk entries', n * n - len(cons))
t = 7
A = np.zeros((len(cons), 2 * t * n))
for row, (a, b) in enumerate(cons):
    for q in range(t):
        A[row, q * n + a] += Phi[q, b]              # d Psi_{q a} * Phi_{q b}
        A[row, t * n + q * n + b] += Psi[q, a]      # Psi_{q a} * d Phi_{q b}
s = np.linalg.svd(A, compute_uv=False)
rk = int((s > 1e-9 * s[0]).sum())
print('variables', 2 * t * n, 'rank', rk, 'tangent dim', 2 * t * n - rk, '| smallest kept sv', s[rk - 1], 'next', s[rk] if rk < len(s) else None)
