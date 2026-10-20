#!/usr/bin/env python3
"""s9_levigroup.py (t2-block3): do the 168 canonical Levi lifts c_A = P_A diag(i^q_A) form a group (exactly, or up to
scalars / Paulis)?  Prints the set of correction factors c_A c_B c_AB^-1."""
import itertools
import numpy as np
from lib3 import *
from collections import Counter

def diag_q(Sm):
    d = []
    for x in range(8):
        v = [(x >> j) & 1 for j in range(3)]
        q = sum(Sm[j][j] * v[j] for j in range(3)) + 2 * sum(Sm[j][k] * v[j] * v[k] for j in range(3) for k in range(j + 1, 3))
        d.append(1j ** (q % 4))
    return np.diag(d)
syms = [[[b[0], b[3], b[4]], [b[3], b[1], b[5]], [b[4], b[5], b[2]]] for b in itertools.product([0, 1], repeat=6)]
def fixes(S):
    if S[:3, 3:].any(): return False
    for v in range(1, 8):
        vec = np.array([(v >> j) & 1 for j in range(3)] * 2)
        w = S @ vec % 2
        if not (w[:3] == w[3:]).all(): return False
    return True
GL = all_gl3()
lift = {}
for A in GL:
    PA = perm_mat(lambda x, A=A: apply_lin(A, x))
    for Sm in syms:
        c = PA @ diag_q(Sm)
        s = sympl(c)
        if s is not None and fixes(s):
            lift[tuple(map(tuple, A))] = c; break
def mul(A, B):
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(3)) % 2 for j in range(3)) for i in range(3))
corr = Counter()
for A in lift:
    for Bm in lift:
        C = mul(A, Bm)
        r = lift[A] @ lift[Bm] @ np.linalg.inv(lift[C])
        # classify: scalar? Pauli times scalar?
        kind = None
        if np.allclose(r, r[0, 0] * I8):
            kind = "scalar %s" % np.round(r[0, 0], 6)
        else:
            for a in range(8):
                for b in range(8):
                    P = pauli(a, b)
                    t = np.trace(P.conj().T @ r) / 8
                    if abs(abs(t) - 1) < 1e-9:
                        kind = "Pauli X%dZ%d x %s" % (a, b, np.round(t, 6))
            if kind is None:
                kind = "other"
        corr[kind] += 1
for k, v in corr.most_common():
    print("count %d: %s" % (v, k))
