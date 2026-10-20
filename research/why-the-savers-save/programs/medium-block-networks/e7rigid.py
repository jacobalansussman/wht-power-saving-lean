#!/usr/bin/env python3
"""cx-medium / e7rigid.py (own code; numpy from our pylib; python3 -I).  Is the rank-7 E7 decoder locally unique?
A 7-total decoder is M = Psi^T Phi (Psi: scatter, Phi: gather, both 7 x 63) with M_aa = 1 and M_ab = 0 on every
commuting pair (anything on the 1,008 anticommuting pairs).  The root decoder M0 = Gram/2 (Phi0 = Psi0 = roots/sqrt2)
is one.  Obvious deformations: GL(7) on (Psi, Phi) (49 dims, M unchanged) and the diagonal rescaling
M -> D^-1 M D (62 dims; it keeps every root support).  If the tangent space of the constraint set at (Psi0, Phi0) has
dimension exactly 49 + 62 = 111, every 7-total decoder near the root decoder is a rescaled root decoder, so its totals
gather sum_b (f.r_b) d_b x_b: supports of at least 27 ports (63 minus the 36 lines of an E6 hyperplane), labels
spanning at least 6 bits, cost >= 6 per total, >= 42 per call = the one-pair break-even.  (Local statement only: other
components of the rank-7 solution set are not excluded.)"""
import sys, os, json
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from designs_cx import design, labels, ip
L = design('E7'); labs, k, B = labels(L); n = len(L)
R = np.array(L, float)                                   # 63 x 8, inside the 7-dim span
# orthonormal basis of the span
U, s, Vt = np.linalg.svd(R, full_matrices=False)
basis = Vt[:7]                                           # 7 x 8
Phi0 = basis @ R.T / np.sqrt(ip(L[0], L[0]) / 2) / np.sqrt(2)   # 7 x 63, so that Phi0^T Phi0 = Gram/2 in norm-2 units
M0 = Phi0.T @ Phi0
cons = [(a, a) for a in range(n)] + [(a, b) for a in range(n) for b in range(n) if a != b and B[a][b] == 1]
assert all(abs(M0[a, b] - (1.0 if a == b else 0.0)) < 1e-9 for a, b in cons), 'root decoder check'
# Jacobian of (Psi, Phi) -> M[cons] at (Phi0, Phi0): dM_ab = sum_i dPsi[i,a] Phi0[i,b] + Psi0[i,a] dPhi[i,b]
nv = 2 * 7 * n
J = np.zeros((len(cons), nv))
for r, (a, b) in enumerate(cons):
    for i in range(7):
        J[r, i * n + a] += Phi0[i, b]                  # dPsi[i, a]
        J[r, 7 * n + i * n + b] += Phi0[i, a]          # dPhi[i, b]  (Psi0 = Phi0)
rk = np.linalg.matrix_rank(J, tol=1e-8)
kern = nv - rk
print(f"constraints {len(cons)}, variables {nv}, Jacobian rank {rk}, tangent dimension {kern}; obvious deformations 49 + 62 = 111")
print('locally rigid (only GL(7) and diagonal rescaling):', kern == 111)
json.dump(dict(constraints=len(cons), variables=nv, rank=int(rk), tangent=int(kern), obvious=111, rigid=bool(kern == 111)),
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'e7rigid.json'), 'w'), indent=1)
