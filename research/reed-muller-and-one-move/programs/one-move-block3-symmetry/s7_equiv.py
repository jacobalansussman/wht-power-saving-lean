#!/usr/bin/env python3
"""s7_equiv.py (t2-block3): checks behind PROOF.md Lemma S6 (only two distinct per-character problems for the Singer).
(a) an order-3 Levi element n with n c n^-1 = c^2 exactly (so n' = K^-1 n K is free and n' c' n'^-1 = c'^2);
(b) complex conjugation: conj(c') = p c' p^-1 for a Pauli p (up to a scalar of order 1), and conj of every engine
    generator is a free letter times (at most) one paid letter; conj(K) = X^7 K.
"""
import itertools
import numpy as np
from lib3 import *

c = np.load("singer_c.npy"); cp = Kinv @ np.linalg.inv(c) @ K
c2 = c @ c
# (a) search Levi lifts P_A diag(i^q) (times Paulis) of order-3 elements A with n c n^-1 = c^2
def diag_q(Sm):
    d = []
    for x in range(8):
        v = [(x >> j) & 1 for j in range(3)]
        q = sum(Sm[j][j] * v[j] for j in range(3)) + 2 * sum(Sm[j][k] * v[j] * v[k] for j in range(3) for k in range(j + 1, 3))
        d.append(1j ** (q % 4))
    return np.diag(d)
syms = [[[b[0], b[3], b[4]], [b[3], b[1], b[5]], [b[4], b[5], b[2]]] for b in itertools.product([0, 1], repeat=6)]
found = None
for A in all_gl3():
    PA = perm_mat(lambda x, A=A: apply_lin(A, x))
    for Sm in syms:
        base = PA @ diag_q(Sm)
        for a in range(8):
            for b in range(8):
                nn = base @ pauli(a, b)
                lhs = nn @ c @ np.linalg.inv(nn)
                # equality up to a scalar
                k0 = np.argmax(np.abs(c2.flatten()))
                lam = lhs.flatten()[k0] / c2.flatten()[k0]
                if np.allclose(lhs, lam * c2) and in_B(K @ nn @ Kinv) and in_B(nn):
                    found = (A, nn, lam); break
            if found: break
        if found: break
    if found: break
A, nn, lam = found
npr = Kinv @ nn @ K
print("(a) n found: A =", A, " n c n^-1 = %s * c^2;" % np.round(lam, 6), "n in B:", in_B(nn), " n' = K^-1 n K in B:", in_B(npr))
lhs = npr @ cp @ np.linalg.inv(npr)
k0 = np.argmax(np.abs((cp @ cp).flatten()))
print("    n' c' n'^-1 = c'^2 up to the scalar", np.round(lhs.flatten()[k0] / (cp @ cp).flatten()[k0], 6),
      ":", np.allclose(lhs, lhs.flatten()[k0] / (cp @ cp).flatten()[k0] * (cp @ cp)))
# (b) complex conjugation
cc = np.conj(cp)
hit = []
for a in range(8):
    for b in range(8):
        p = pauli(a, b)
        t = p @ cp @ np.linalg.inv(p)
        k0 = np.argmax(np.abs(cc.flatten()))
        sc = cc.flatten()[k0] / t.flatten()[k0]
        if np.allclose(cc, sc * t):
            hit.append(((a, b), np.round(sc, 6)))
print("(b) conj(c') = s p c' p^-1 for Paulis p = X^a Z^b:", hit[:4], "... (%d found)" % len(hit))
print("    conj(K) == X^7 K:", np.allclose(np.conj(K), X(7) @ K))
for z in range(1, 8):
    assert np.allclose(np.conj(dirz(z)), X(z) @ dirz(z))
print("    conj(dir z) == X^z dir z for z = 1..7: True")
for g in [X(1), Z(3), phase_letter(5), phase_letter(2, 1), perm_mat(lambda x: x ^ ((x & 1) << 2))]:
    assert in_B(np.conj(g))
print("    conj of the free generators is free: True")
