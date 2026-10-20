#!/usr/bin/env python3
"""s2_split.py (t2-block3): eigen-splits that decide the per-character hook test of symmetric block-3 shapes.

For a symmetry element c of exact order n acting through the Levi (c and c' = K^-1 c^-1 K in B):
  - demand split: multiplicities of the n-th roots of unity of rho(c') on N/E_N (a fixed array's projected pieces);
  - cap split of a FIXED move m (c m c^-1 = m): multiplicities of rho(c) on (N + rho(m) N)/N;
  - fixed cells / fixed moves of c in the dual polar graph (Lagrangians of F_2^6), from symplectic images.
Traces are computed from matrix powers c^s directly (no accumulated products).
"""
import itertools
import numpy as np
from lib3 import *

out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

Nb = np.zeros((560, 280), dtype=complex); Nb[NIDX, np.arange(280)] = 1
RKi = rho(Kinv)
EN = intersect(Nb, colspace(RKi @ Nb))
assert EN.shape[1] == 56

def tr_on(basis, R):
    return np.trace(np.linalg.pinv(basis) @ R @ basis)

def mults(trs, n):
    z = np.exp(2j * np.pi / n)
    res = []
    for k in range(n):
        m = sum(trs[s] * z ** (-k * s) for s in range(n)) / n
        assert abs(m.imag) < 1e-6 and abs(m.real - round(m.real)) < 1e-6, m
        res.append(int(round(m.real)))
    return res

def demand_split(cp, n):
    trs = []
    for s in range(n):
        R = rho(np.linalg.matrix_power(cp, s))
        trs.append(tr_on(Nb, R) - tr_on(EN, R))
    return mults(trs, n)

def cap_split(c, m, n):
    assert np.allclose(c @ m, m @ c), "move does not commute with c"
    Rm = rho(m)
    NM = colspace(np.hstack([Nb, Rm @ Nb]))
    trs = []
    for s in range(n):
        R = rho(np.linalg.matrix_power(c, s))
        trs.append(tr_on(NM, R) - tr_on(Nb, R))
    return NM.shape[1] - 280, mults(trs, n)

# ---------------------------------------------------------------- Lagrangians of F_2^6 and the dual polar graph
def form(u, v):  # symplectic form on (x|z) labels as 6-bit ints: bits 0-2 = x, bits 3-5 = z
    xu, zu, xv, zv = u & 7, u >> 3, v & 7, v >> 3
    return (bin(xu & zv).count("1") + bin(zu & xv).count("1")) & 1

def span(vs):
    S = {0}
    for v in vs:
        S |= {s ^ v for s in S}
    return frozenset(S)

lags = set()
for a, b, c3 in itertools.combinations(range(1, 64), 3):
    if form(a, b) or form(a, c3) or form(b, c3):
        continue
    S = span([a, b, c3])
    if len(S) == 8:
        lags.add(S)
lags = sorted(lags, key=lambda S: sorted(S))
say("Lagrangians:", len(lags))
LZ = span([8, 16, 32])
idx = {S: i for i, S in enumerate(lags)}
edges = [(i, j) for i in range(len(lags)) for j in range(i + 1, len(lags)) if len(lags[i] & lags[j]) == 4]
say("dual polar graph edges:", len(edges))

def act(Sm, S):
    """Sm: 6x6 F2 matrix on coordinates (x0,x1,x2,z0,z1,z2)"""
    res = set()
    for v in S:
        vec = np.array([(v >> j) & 1 for j in range(6)])
        w = Sm @ vec % 2
        res.add(int(sum(int(w[j]) << j for j in range(6))))
    return frozenset(res)

def fixed_stats(g, name):
    Sm = sympl(g)
    fc = [S for S in lags if act(Sm, S) == S]
    fe = [(i, j) for (i, j) in edges if act(Sm, lags[i]) == lags[i] and act(Sm, lags[j]) == lags[j]]
    fe_sw = [(i, j) for (i, j) in edges if act(Sm, lags[i]) == lags[j] and act(Sm, lags[j]) == lags[i]]
    dist = lambda S: 3 - (len(S & LZ).bit_length() - 1)
    say("%s: fixed cells %d (distances from L_Z: %s), fixed edges %d, swapped edges %d" % (
        name, len(fc), sorted(dist(S) for S in fc), len(fe), len(fe_sw)))

# ---------------------------------------------------------------- symmetry elements
def coordperm(p):   # x -> bits permuted: bit j of x goes to bit p[j]
    def f(x):
        y = 0
        for j in range(3):
            if (x >> j) & 1:
                y |= 1 << p[j]
        return y
    return perm_mat(f)

c3 = coordperm([1, 2, 0])
assert np.allclose(K @ c3 @ Kinv, c3)
c2 = coordperm([1, 0, 2])
assert np.allclose(K @ c2 @ Kinv, c2)
c7 = np.load("singer_c.npy")
say("coordinate permutations commute with K: ok")

fixed_stats(c7, "Singer c (order 7)")
fixed_stats(c3, "coordinate 3-cycle (order 3)")
fixed_stats(c2, "coordinate transposition (order 2)")
for v in (1, 3, 7):
    fixed_stats(X(v), "Pauli X^%d (order 2)" % v)

# ---------------------------------------------------------------- demand splits
for name, c, n in [("Singer", c7, 7), ("3-cycle", c3, 3), ("transposition", c2, 2), ("X^1", X(1), 2), ("X^7", X(7), 2)]:
    cp = Kinv @ np.linalg.inv(c) @ K
    assert in_B(cp)
    say("demand split of c' on N/E_N (%s, order %d): %s" % (name, n, demand_split(cp, n)))

# ---------------------------------------------------------------- cap splits of fixed moves
S_all = np.diag([1j ** (bin(x).count("1")) for x in range(8)])        # i^{wt x}: in B, commutes with coordinate perms
CZ_all = np.diag([(-1) ** ((x & 1) * ((x >> 1) & 1) + ((x >> 1) & 1) * ((x >> 2) & 1) + (x & 1) * ((x >> 2) & 1)) for x in range(8)]).astype(complex)
moves3 = {"dir7": dirz(7), "S dir7 S^-1": S_all @ dirz(7) @ np.linalg.inv(S_all),
          "CZ dir7 CZ": CZ_all @ dirz(7) @ CZ_all, "Z^7 dir7 Z^7": Z(7) @ dirz(7) @ Z(7)}
for nm, m in moves3.items():
    say("3-cycle fixed move %-14s cap, split over chars of c: %s" % (nm, cap_split(c3, m, 3)))
for z in range(1, 8):
    say("X^1 fixed move dir%d: cap, split: %s" % (z, cap_split(X(1), dirz(z), 2)))
for nm, m in [("dir7", dirz(7)), ("dir3", dirz(3)), ("dir4", dirz(4))]:
    if np.allclose(c2 @ m, m @ c2):
        say("transposition fixed move %s: cap, split: %s" % (nm, cap_split(c2, m, 2)))
open("out_s2_split.txt", "w").write("\n".join(out) + "\n")
