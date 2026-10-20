#!/usr/bin/env python3
"""s6_rf.py (t2-block3): the filtered rung (ladder Theorem RF) in the f = 5 hook module, for full arrays and for
Singer-projected pieces, to screen the per-character programs of the Z7 survivors.

V_0 = N, V_(j+1) = B-span(V_j + dir1 V_j) (keeps every free letter, D V_j inside V_(j+1)).  For a per-character program
with A full arrays and P pieces of character chi:  T >= j + (A dem_j + P dem_j(chi)) / 80  for every j with dem_j > 0
(Theorem RF applied with inputs restricted to N (+) ... (+) N_chi, N_chi the chi-eigenspace of c' on N; the proof of
RF only uses the potential on the input subspace, so restricting inputs to a c'-eigenspace changes only the demands).
Calibration: for standard units (P = 0) the ladder says RF adds nothing at 5 in 14 / 6 in 17.
"""
import numpy as np
from lib3 import *

def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
out = []
Nb = np.zeros((560, 280), dtype=complex); Nb[NIDX, np.arange(280)] = 1
gens = [X(1), X(2), X(4), Z(1), Z(2), Z(4), phase_letter(1), phase_letter(2), phase_letter(4)]
for i in range(3):
    for j in range(3):
        if i != j:
            gens.append(perm_mat(lambda x, i=i, j=j: x ^ (((x >> i) & 1) << j)))
Rg = [rho(g) for g in gens]
for g in gens:
    assert in_B(g)

def bspan(U):
    while True:
        V = colspace(np.hstack([U] + [R @ U for R in Rg]), tol=1e-8)
        if V.shape[1] == U.shape[1]:
            return V
        U = V

RD = rho(dirz(1)); RK = rho(K)
V = [colspace(Nb)]
for j in range(3):
    V.append(bspan(colspace(np.hstack([V[-1], RD @ V[-1]]), tol=1e-8)))
say("dims V_0..V_3:", [v.shape[1] for v in V])
f = colspace(np.hstack([Nb, RD @ Nb])).shape[1] - 280
say("cap f =", f)

c = np.load("singer_c.npy"); cp = Kinv @ np.linalg.inv(c) @ K
Rc = rho(cp)
# chi-eigenspaces of c' on N
z7 = np.exp(2j * np.pi / 7)
Nchi = []
for k in range(7):
    P = sum(z7 ** (-k * h) * np.linalg.matrix_power(Rc, h) for h in range(7)) / 7
    Nchi.append(colspace(P @ Nb, tol=1e-8))
say("dims of the c'-eigenspaces of N:", [x.shape[1] for x in Nchi])

def dem(U, j):
    """dim U - dim { v in U : K v in V_j }"""
    Vj = V[j]
    # K v in V_j  <=>  (I - Vj Vj^*) K U a = 0
    Mx = (np.eye(560) - Vj @ Vj.conj().T) @ RK @ U
    s = np.linalg.svd(Mx, compute_uv=False)
    r = int((s > 1e-8 * max(1, s[0])).sum())
    return r

for j in range(4):
    say("j = %d: dem_j(full array) = %d; dem_j(piece) per character = %s" % (j, dem(Nb, j), [dem(Nchi[k], j) for k in range(7)]))
open("out_s6_rf.txt", "w").write("\n".join(out) + "\n")
