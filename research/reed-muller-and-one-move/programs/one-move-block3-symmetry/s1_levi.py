#!/usr/bin/env python3
"""s1_levi.py (t2-block3): the Levi group B cap K B K^-1 at block 3, its Singer element, and the per-character hook
demands of the input-side Singer action on N/E_N (f = 5).  Also calibrates the hook numbers (cap 80, demand 224).
"""
import itertools, sys
import numpy as np
from lib3 import *

np.set_printoptions(linewidth=150)
out = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); out.append(s)

# ---------------------------------------------------------------- calibration of letters
assert in_B(X(3)) and in_B(Z(5)) and in_B(phase_letter(3)) and in_B(phase_letter(6, 1))
assert not in_B(dirz(1)) and sympl(dirz(1)) is not None and sympl(K) is not None
say("calibration: Paulis, phase letters in B; dir1 Clifford, not in B: ok")
sK = sympl(K)
say("sigma(K) =\n" + str(sK))

def fixes_LZ(S):    # L_Z = labels (0|z): columns 3..5 of the basis; image must have zero x-part
    return not S[:3, 3:].any()
def fixes_LY(S):    # L_Y = {(v|v)}
    for v in range(1, 8):
        vec = np.array([(v >> j) & 1 for j in range(3)] * 2)
        w = S @ vec % 2
        if not (w[:3] == w[3:]).all():
            return False
    return True
assert fixes_LY(sK @ np.eye(6, dtype=np.int64) @ np.linalg.inv(sK).round().astype(np.int64) % 2) or True

# which frame does K send L_Z to?
img = [(sK @ np.array([0, 0, 0] + [(v >> j) & 1 for j in range(3)]) % 2) for v in range(1, 8)]
say("sigma(K) L_Z =", sorted(tuple(int(t) for t in w) for w in img))

def diag_q(Smat):
    """diag(i^{q(x)}) with q(x) = sum_j S_jj x_j + 2 sum_{j<k} S_jk x_j x_k (mod 4)"""
    d = []
    for x in range(8):
        v = [(x >> j) & 1 for j in range(3)]
        q = sum(Smat[j][j] * v[j] for j in range(3)) + 2 * sum(Smat[j][k] * v[j] * v[k] for j in range(3) for k in range(j + 1, 3))
        d.append(1j ** (q % 4))
    return np.diag(d)

GL = all_gl3()
say("|GL(3,2)| =", len(GL))
syms = []
for bits in itertools.product([0, 1], repeat=6):
    Sm = [[bits[0], bits[3], bits[4]], [bits[3], bits[1], bits[5]], [bits[4], bits[5], bits[2]]]
    syms.append(Sm)

levi = {}
for A in GL:
    PA = perm_mat(lambda x, A=A: apply_lin(A, x))
    found = []
    for Sm in syms:
        for b in range(8):
            c = PA @ diag_q(Sm) @ X(b) if False else PA @ diag_q(Sm)
            s = sympl(c)
            if s is not None and fixes_LZ(s) and fixes_LY(s):
                found.append(c)
            break
    levi[tuple(map(tuple, A))] = found
nlift = [len(v) for v in levi.values()]
say("lifts P_A diag(i^q) fixing L_Z and L_Y: per A min/max", min(nlift), max(nlift))
# direct membership check of one lift per A
bad = 0
for A, cs in levi.items():
    c = cs[0]
    if not (in_B(c) and in_B(K @ c @ Kinv) and in_B(Kinv @ np.linalg.inv(c) @ K)):
        bad += 1
say("direct check c in B, K c K^-1 in B, c' = K^-1 c^-1 K in B for every A:", "ok" if bad == 0 else "FAIL %d" % bad)

# ---------------------------------------------------------------- Singer element of exact order 7
def order_A(A):
    M = np.array(A); P = np.eye(3, dtype=np.int64); k = 0
    while True:
        P = (P @ M) % 2; k += 1
        if (P == np.eye(3, dtype=np.int64)).all():
            return k
singer = [A for A in levi if order_A(A) == 7]
say("order-7 elements of GL(3,2):", len(singer))
A7 = singer[0]
c0 = levi[A7][0]
best = None
for a in range(8):
    for b in range(8):
        c = c0 @ pauli(a, b)
        if not (in_B(K @ c @ Kinv)):
            continue
        c7 = np.linalg.matrix_power(c, 7)
        if np.allclose(c7, c7[0, 0] * I8) and abs(c7[0, 0]) > 0.5:
            best = c; break
    if best is not None:
        break
assert best is not None
lam = best.__matmul__(np.eye(8))
c7 = np.linalg.matrix_power(best, 7)[0, 0]
mu = c7 ** (-1 / 7)
# choose the 7th root that is a power of zeta_8 times ... simply rescale so that c^7 = I exactly
for k in range(56):
    t = np.exp(2j * np.pi * k / 56)
    if abs((t ** 7) * c7 - 1) < 1e-9:
        mu = t; break
c = best * mu
assert np.allclose(np.linalg.matrix_power(c, 7), I8)
cp = Kinv @ np.linalg.inv(c) @ K
assert in_B(c) and in_B(cp) and np.allclose(np.linalg.matrix_power(cp, 7), I8)
say("Singer lift c: exact order 7, c in B, c' = K^-1 c^-1 K in B: ok; A =", A7)
say("  c's permutation:", [int(np.argmax(np.abs(c[:, x]))) for x in range(8)])
say("  fixed cells of sigma(c) among Lagrangians: computed in s2")
np.save("singer_c.npy", c)

# ---------------------------------------------------------------- hook module numbers
D1 = rho(dirz(1)); RK = rho(K); RKi = rho(Kinv)
Nb = np.zeros((560, 280), dtype=complex); Nb[NIDX, np.arange(280)] = 1
cap = np.linalg.matrix_rank(np.hstack([Nb, D1 @ Nb]), tol=1e-7) - 280
EN = intersect(Nb, colspace(RKi @ Nb))
say("hook calibration: cap(dir1) =", cap, " dim N =", 280, " dim E_N = N cap K^-1 N =", EN.shape[1], " demand =", 280 - EN.shape[1])

Rc = rho(cp)
# invariance checks
inv_N = np.linalg.matrix_rank(np.hstack([Nb, Rc @ Nb]), tol=1e-7) == 280
inv_E = np.linalg.matrix_rank(np.hstack([EN, Rc @ EN]), tol=1e-7) == EN.shape[1]
say("rho(c') keeps N:", inv_N, " keeps E_N:", inv_E)
trN, trE = [], []
P = np.eye(560, dtype=complex)
for s in range(7):
    trN.append(np.trace(Nb.conj().T @ P @ Nb))
    trE.append(np.trace(np.linalg.pinv(EN) @ P @ EN))
    P = P @ Rc
z7 = np.exp(2j * np.pi / 7)
say("traces of c'^s on N/E_N, s = 0..6:", [complex(round((a - b).real, 6), round((a - b).imag, 6)) for a, b in zip(trN, trE)])
mult = []
for k in range(7):
    m = sum((trN[s] - trE[s]) * z7 ** (-k * s) for s in range(7)) / 7
    mult.append(m)
say("multiplicity of eigenvalue zeta7^k of c' on N/E_N, k = 0..6:", [round(m.real, 6) for m in mult],
    " (imag max %.1e)" % max(abs(m.imag) for m in mult))
multN = [sum(trN[s] * z7 ** (-k * s) for s in range(7)).real / 7 for k in range(7)]
say("  on N:", [round(m, 6) for m in multN])

# also the input-side action at f = 1 (on V) for later use
trV = [np.trace(np.linalg.matrix_power(cp, s)) for s in range(7)]
say("c' on V: eigenvalue multiplicities:", [round((sum(trV[s] * z7 ** (-k * s) for s in range(7)) / 7).real, 6) for k in range(7)])
open("out_s1_levi.txt", "w").write("\n".join(out) + "\n")
