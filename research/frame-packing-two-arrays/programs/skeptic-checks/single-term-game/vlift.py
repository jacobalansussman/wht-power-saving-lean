#!/usr/bin/env python3
"""vlift.py: does an apartment swap (an Sp(2,2) element exchanging two of the three lines Z, X, Y) have an
involutive lift in the engine's one-qubit group G = <Paulis, S, scalars i^k, C>?  Lifts of a given Sp-image are
n*g0 with n in N = Paulis x {1, i, -1, -i}; we list (n g0)^2 for all 16 n."""
import numpy as np, itertools
I = np.eye(2, dtype=complex); X = np.array([[0, 1], [1, 0]], dtype=complex); Z = np.diag([1, -1]).astype(complex)
Y = 1j * X @ Z; S = np.diag([1, 1j]); C = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
H = np.array([[1, 1], [1, -1]]) / np.sqrt(2); h = np.exp(1j * np.pi / 4) * H
N = [c * P for c in (1, 1j, -1, -1j) for P in (I, X, Y, Z)]
def in_Qi(M):   # entries of the form (a + b i)/2^k with small rational a, b? test: real and imag parts times 2^4 integral
    return np.allclose(np.round(M * 16), M * 16)
for name, g0 in (("Z<->X (h = e^{i pi/4} H)", h), ("Z<->Y (C)", C), ("X<->Y (S)", S)):
    sq = [n @ g0 @ n @ g0 for n in N]
    inv = [np.allclose(q, I) for q in sq]
    print(name, "g0 entries in Q(i):", in_Qi(g0), "| any involutive lift:", any(inv),
          "| squares are scalar multiples of I:", sum(np.allclose(q, q[0, 0] * I) for q in sq), "of 16")
print("H entries in Q(i):", in_Qi(H))
