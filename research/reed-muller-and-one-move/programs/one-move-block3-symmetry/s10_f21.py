#!/usr/bin/env python3
"""s10_f21.py (t2-block3): the F21 refinement of the trivial-character problem of the Z7 survivor 15 in 42.

F21 = <c, n>, c the Singer lift (order 7), n an order-3 Levi element with n c n^-1 = c^2 (s7_equiv.py).  Array side:
c'_g = K^-1 g^-1 K.  In the Z7-trivial problem the residual Z3 = <n> fixes all three inputs (two orbit sums and the
trivial piece), so they split by characters of <n>; moves of F21-orbits of size 7 are <n>-fixed (projected) moves.
Computes: d = split of rho(c'_n) on N/E_N; p = split on the c'-invariant part of N/E_N (dim 32); sigma = split of
rho(n) on (N + rho(m)N)/N for move letters m commuting with n.  Then solves the tight equations of PROOF.md S9.
"""
import itertools
import numpy as np
from lib3 import *

out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

c = np.load("singer_c.npy")
# rebuild n as in s7_equiv.py
def diag_q(Sm):
    d = []
    for x in range(8):
        v = [(x >> j) & 1 for j in range(3)]
        q = sum(Sm[j][j] * v[j] for j in range(3)) + 2 * sum(Sm[j][k] * v[j] * v[k] for j in range(3) for k in range(j + 1, 3))
        d.append(1j ** (q % 4))
    return np.diag(d)
syms = [[[b[0], b[3], b[4]], [b[3], b[1], b[5]], [b[4], b[5], b[2]]] for b in itertools.product([0, 1], repeat=6)]
c2 = c @ c
n = None
for A in all_gl3():
    PA = perm_mat(lambda x, A=A: apply_lin(A, x))
    for Sm in syms:
        base = PA @ diag_q(Sm)
        for a in range(8):
            for b in range(8):
                g = base @ pauli(a, b)
                if np.allclose(g @ c @ np.linalg.inv(g), c2) and in_B(K @ g @ Kinv) and in_B(g):
                    n = g; break
            if n is not None: break
        if n is not None: break
    if n is not None: break
n3 = np.linalg.matrix_power(n, 3)
lam = n3[0, 0]
assert np.allclose(n3, lam * I8)
for k in range(24):
    t = np.exp(2j * np.pi * k / 24)
    if abs(t ** 3 * lam - 1) < 1e-9:
        n = n * t; break
assert np.allclose(np.linalg.matrix_power(n, 3), I8) and np.allclose(n @ c @ np.linalg.inv(n), c2)
say("n: exact order 3, n c n^-1 = c^2, n and K n K^-1 free")
cp = Kinv @ np.linalg.inv(c) @ K
cpn = Kinv @ np.linalg.inv(n) @ K
assert in_B(cpn)

Nb = np.zeros((560, 280), dtype=complex); Nb[NIDX, np.arange(280)] = 1
EN = intersect(Nb, colspace(rho(Kinv) @ Nb)); assert EN.shape[1] == 56
z3 = np.exp(2j * np.pi / 3)

def mults(trs, nn):
    z = np.exp(2j * np.pi / nn); res = []
    for k in range(nn):
        m = sum(trs[s] * z ** (-k * s) for s in range(nn)) / nn
        assert abs(m.imag) < 1e-6 and abs(m.real - round(m.real)) < 1e-6, m
        res.append(int(round(m.real)))
    return res

R = [rho(np.linalg.matrix_power(cpn, s)) for s in range(3)]
d = mults([np.trace(np.linalg.pinv(Nb) @ Rs @ Nb) - np.trace(np.linalg.pinv(EN) @ Rs @ EN) for Rs in R], 3)
say("d  = split of c'_n on N/E_N:", d)
# c'-invariant parts
Rc = rho(cp)
P0 = sum(np.linalg.matrix_power(Rc, h) for h in range(7)) / 7
N0 = colspace(P0 @ Nb); E0 = colspace(P0 @ EN)
say("dims of c'-invariant parts: N", N0.shape[1], " E_N", E0.shape[1])
p = mults([np.trace(np.linalg.pinv(N0) @ Rs @ N0) - np.trace(np.linalg.pinv(E0) @ Rs @ E0) for Rs in R], 3)
say("p  = split of c'_n on the c'-invariant part of N/E_N:", p)

# moves commuting with n: search b dir z b^-1 over a menu of free letters b
menu = [I8]
gens = [X(1), X(2), X(4), Z(1), Z(2), Z(4), phase_letter(1), phase_letter(2), phase_letter(4), phase_letter(7)]
for g in gens:
    menu.append(g)
for g, h in itertools.product(gens, gens):
    menu.append(g @ h)
Rn = [rho(np.linalg.matrix_power(n, s)) for s in range(3)]
found = {}
for z in range(1, 8):
    for b in menu:
        m = b @ dirz(z) @ np.linalg.inv(b)
        r = n @ m @ np.linalg.inv(n)
        k0 = np.argmax(np.abs(m.flatten()))
        sc = r.flatten()[k0] / m.flatten()[k0]
        if np.allclose(r, sc * m):
            NM = colspace(np.hstack([Nb, rho(m) @ Nb]))
            sg = mults([np.trace(np.linalg.pinv(NM) @ Rs @ NM) - np.trace(Nb.conj().T @ Rs @ Nb) for Rs in Rn], 3) if np.isclose(sc, 1) else None
            found.setdefault((z, np.round(sc, 6)), sg)
            break
for (z, sc), sg in found.items():
    say("move b dir%d b^-1 commuting with n (n m n^-1 = %s m): sigma = %s" % (z, sc, sg))

# tight equations for the Z7-trivial problem of 15 in 42 under <n>: r free <n>-orbits (F21-orbits of 21 moves),
# q <n>-fixed moves (F21-orbits of 7), 3r + q = 6; per <n>-character k (array side eigenvalue z3^k):
#   80 r + q * sigma_{-k}  ==  2 d_k + p_k
sols = []
for sg in set(tuple(v) for v in found.values() if v is not None):
    for r_ in range(3):
        q_ = 6 - 3 * r_
        if all(80 * r_ + q_ * sg[(-k) % 3] == 2 * d[k] + p[k] for k in range(3)):
            sols.append((sg, r_, q_))
say("tight solutions (sigma, r, q):", sols if sols else "NONE")
open("out_s10_f21.txt", "w").write("\n".join(out) + "\n")
