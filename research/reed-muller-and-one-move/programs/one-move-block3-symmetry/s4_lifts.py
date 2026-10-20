#!/usr/bin/env python3
"""s4_lifts.py (t2-block3): demand and fixed-move cap splits for EVERY Pauli-twisted lift c P of the 3-cycle and of the
transposition (P runs over the 64 Paulis; c P is rescaled to exact order n when (c P)^n is a scalar), and for the
Singer lift's Pauli twists of exact order 7.  Fixed moves tested: dir z (z = 1..7) and the B-conjugates S dir z S^-1,
(S = i^{wt x}), whenever they commute with the lift.
"""
import numpy as np
from lib3 import *
import itertools
Nb = np.zeros((560, 280), dtype=complex); Nb[NIDX, np.arange(280)] = 1
EN = intersect(Nb, colspace(rho(Kinv) @ Nb)); assert EN.shape[1] == 56
pN, pE = np.linalg.pinv(Nb), np.linalg.pinv(EN)
def mults(trs, n):
    z = np.exp(2j * np.pi / n); res = []
    for k in range(n):
        m = sum(trs[s] * z ** (-k * s) for s in range(n)) / n
        assert abs(m.imag) < 1e-6 and abs(m.real - round(m.real)) < 1e-6, m
        res.append(int(round(m.real)))
    return res
def coordperm(p):
    def f(x):
        y = 0
        for j in range(3):
            if (x >> j) & 1: y |= 1 << p[j]
        return y
    return perm_mat(f)
NMc = {}
def NMof(mn, m):
    if mn not in NMc:
        NM = colspace(np.hstack([Nb, rho(m) @ Nb])); NMc[mn] = (NM, np.linalg.pinv(NM))
    return NMc[mn]

c3 = coordperm([1, 2, 0]); c2 = coordperm([1, 0, 2]); c7 = np.load("singer_c.npy")
S_all = np.diag([1j ** (bin(x).count("1")) for x in range(8)])
CZ_all = np.diag([(-1) ** ((x & 1) * ((x >> 1) & 1) + ((x >> 1) & 1) * ((x >> 2) & 1) + (x & 1) * ((x >> 2) & 1)) for x in range(8)]).astype(complex)
cand_moves = {}
for z in range(1, 8):
    cand_moves["dir%d" % z] = dirz(z)
    cand_moves["S dir%d S'" % z] = S_all @ dirz(z) @ np.linalg.inv(S_all)
lines = []
for name, c, n in [("3-cycle", c3, 3), ("transposition", c2, 2), ("Singer", c7, 7)]:
    seen = {}
    for a in range(8):
        for b in range(8):
            g = c @ pauli(a, b)
            gn = np.linalg.matrix_power(g, n)
            if not np.allclose(gn, gn[0, 0] * I8):
                continue
            lam = gn[0, 0]
            for k in range(8 * n):
                t = np.exp(2j * np.pi * k / (8 * n))
                if abs(t ** n * lam - 1) < 1e-9:
                    g = g * t; break
            assert np.allclose(np.linalg.matrix_power(g, n), I8)
            gp = Kinv @ np.linalg.inv(g) @ K
            assert in_B(g) and in_B(gp)
            Rg = [rho(np.linalg.matrix_power(g, k)) for k in range(n)]
            Rgp = [rho(np.linalg.matrix_power(gp, k)) for k in range(n)]
            dsp = tuple(mults([np.trace(pN @ R @ Nb) - np.trace(pE @ R @ EN) for R in Rgp], n))
            caps = []
            for mn, m in cand_moves.items():
                if np.allclose(g @ m, m @ g):
                    NM, pNM = NMof(mn, m)
                    caps.append((mn, tuple(mults([np.trace(pNM @ R @ NM) - np.trace(pN @ R @ Nb) for R in Rg], n))))
            capset = sorted(set(cs for _, cs in caps))
            key = (dsp, tuple(capset))
            seen.setdefault(key, []).append("X%dZ%d" % (a, b))
    for (dsp, capset), lst in seen.items():
        s = "%s lifts (%d of them, e.g. %s): demand split %s; fixed-move cap splits %s" % (name, len(lst), ",".join(lst[:4]), dsp, list(capset) if capset else "none among the tested moves")
        print(s, flush=True); lines.append(s)
open("out_s4_lifts.txt", "w").write("\n".join(lines) + "\n")
