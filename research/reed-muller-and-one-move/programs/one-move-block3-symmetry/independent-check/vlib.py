#!/usr/bin/env python3
"""vlib.py (v-t2-block3, skeptic's own code, written without importing t2-block3's lib3.py).
Block-3 engine letters as 8x8 complex matrices; symplectic images on labels (x|z); Lagrangians of F2^6;
the f = 5 hook module V (x) Lambda^4 V with Y = N = span{e_x (x) e_T : x in T}."""
import itertools, numpy as np
n8 = 8
def Xm(a):
    M = np.zeros((8, 8), complex)
    for x in range(8): M[x ^ a, x] = 1
    return M
def Zm(b): return np.diag([(-1.0) ** bin(b & x).count("1") for x in range(8)]).astype(complex)
def pau(a, b): return Xm(a) @ Zm(b)
ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
def dirm(z): return ca * np.eye(8) + cb * Xm(z)
Kmat = dirm(1) @ dirm(2) @ dirm(4)
Kinv = np.linalg.inv(Kmat)
inv = np.linalg.inv
def lin(A, x):  # A as tuple of 3 column images (ints): x = sum x_j e_j -> sum x_j A[j]
    y = 0
    for j in range(3):
        if (x >> j) & 1: y ^= A[j]
    return y
def permm(f):
    M = np.zeros((8, 8), complex)
    for x in range(8): M[f(x), x] = 1
    return M
def gl3():
    out = []
    for cols in itertools.product(range(1, 8), repeat=3):
        if len({lin(cols, x) for x in range(8)}) == 8: out.append(cols)
    return out
# Pauli label vector: integer 0..63 = a | (b << 3) for X^a Z^b
PAULIS = {}
for a in range(8):
    for b in range(8): PAULIS[(a, b)] = pau(a, b)
def sym(g):
    """6x6 F2 matrix: column j = label of g P_j g^-1 (P_j = X^{e_j}, j<3; Z^{e_{j-3}}, j>=3); None if not Clifford"""
    gi = inv(g); cols = []
    for j in range(6):
        P = PAULIS[(1 << j, 0)] if j < 3 else PAULIS[(0, 1 << (j - 3))]
        Q = g @ P @ gi; hit = None
        for (a, b), R in PAULIS.items():
            t = np.vdot(R, Q) / 8
            if abs(abs(t) - 1) < 1e-8: hit = (a, b); break
        if hit is None: return None
        a, b = hit
        cols.append([(a >> i) & 1 for i in range(3)] + [(b >> i) & 1 for i in range(3)])
    return np.array(cols, dtype=np.int64).T % 2
def is_free(g):
    nz = np.abs(g) > 1e-9
    if not ((nz.sum(0) == 1).all() and (nz.sum(1) == 1).all()): return False
    v = g[nz]; r = v / v[0]
    if not np.allclose(np.abs(r), 1): return False
    q = np.angle(r) / (np.pi / 2)
    if not np.allclose(q, np.round(q), atol=1e-8): return False
    return sym(g) is not None
def vec2int(v): return int(sum(int(v[i]) << i for i in range(6)))
def int2vec(u): return np.array([(u >> i) & 1 for i in range(6)], dtype=np.int64)
def form(u, w):  # symplectic form on labels (x|z): x.z' + z.x'
    ux, uz, wx, wz = u & 7, u >> 3, w & 7, w >> 3
    return (bin(ux & wz).count("1") + bin(uz & wx).count("1")) & 1
def span(vs):
    s = {0}
    for v in vs: s = s | {x ^ v for x in s}
    return frozenset(s)
def lagrangians():
    L = set()
    for tri in itertools.combinations(range(1, 64), 3):
        if all(form(u, w) == 0 for u, w in itertools.combinations(tri, 2)):
            S = span(tri)
            if len(S) == 8: L.add(S)
    return sorted(L, key=lambda s: sorted(s))
def apply_sym(S, sub):
    return frozenset(vec2int(S @ int2vec(u) % 2) for u in sub)
LZ = frozenset(b << 3 for b in range(8))
# ---------------- hook module
S4 = list(itertools.combinations(range(8), 4))
IDX = {T: i for i, T in enumerate(S4)}
def lam4(g):
    U = np.array(S4)
    return np.linalg.det(g[U[:, None, :, None], U[None, :, None, :]])
def rho(g): return np.kron(g, lam4(g))
NCOLS = [x * 70 + IDX[T] for x in range(8) for T in S4 if x in T]
Nb = np.zeros((560, 280), complex); Nb[NCOLS, np.arange(280)] = 1
def orth(A, tol=1e-8):
    if A.shape[1] == 0: return A
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    r = int((s > tol * max(1, s[0])).sum()); return U[:, :r]
def cap_meet(A, B):
    A, B = orth(A), orth(B)
    M = np.hstack([A, -B]); U, s, Vh = np.linalg.svd(M)
    r = int((s > 1e-8).sum()); ns = Vh[r:].conj().T
    return orth(A @ ns[:A.shape[1]])
def quot_trace(Rg, big, small):
    """trace of Rg on big/small (both Rg-stable, orthonormal bases, small inside big)"""
    big, small = orth(big), orth(small)
    return np.trace(big.conj().T @ Rg @ big) - (np.trace(small.conj().T @ Rg @ small) if small.shape[1] else 0)
def mults(trs, order):
    z = np.exp(2j * np.pi / order); out = []
    for k in range(order):
        m = sum(trs[s] * z ** (-k * s) for s in range(order)) / order
        assert abs(m.imag) < 1e-5 and abs(m.real - round(m.real)) < 1e-5, m
        out.append(int(round(m.real)))
    return out
