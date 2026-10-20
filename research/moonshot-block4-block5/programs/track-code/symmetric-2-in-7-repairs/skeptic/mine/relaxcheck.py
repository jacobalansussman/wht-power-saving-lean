"""relaxcheck.py (v-mw2-absorb4): own frame-sum relaxation (exact E frames: cells g, g' in one frame iff g' g^-1 is an
E-free letter) of a dumped program, and the anatomy of array 0's kernel-frame block (relative letter u, its order,
freeness, centrality)."""
import sys, pickle
import numpy as np
from mycheck import *
from scipy.optimize import minimize
dmp = pickle.load(open(sys.argv[1], 'rb')); m = dmp['m']; W = dmp['W']
isfree = make_free_test('exact', m)
gens = dmp['gens']; K = dmp['K']
reps = []   # (l, matrix) frame representatives
def frame_of(l, M):
    for k, (l2, R) in enumerate(reps):
        if l2 == l and isfree(M @ np.conj(R).T):
            return k
    reps.append((l, M)); return len(reps) - 1
kf = frame_of(0, K)
cellframe = [[(None if isfree(h) else frame_of(l, h)) for (l, h) in T] for T in gens]
nfr = len(reps)
print("frames:", nfr, "kernel frame index", kf)
for i, T in enumerate(gens):
    print(f"  g{i}: " + "  ".join(f"x{l}:F{cellframe[i][j]}" for j, (l, h) in enumerate(T)))
def relax(alpha):
    cv = coefs(dmp, alpha)
    A = np.zeros((nfr, len(gens)), dtype=complex)
    for i in range(len(gens)):
        for j in range(len(gens[i])):
            if cellframe[i][j] is not None:
                A[cellframe[i][j], i] += cv[i][j]
    t = np.zeros(nfr, dtype=complex); t[kf] = 1
    b = np.linalg.lstsq(A, t, rcond=None)[0]
    return float(np.sum(np.abs(A @ b - t) ** 2))
al = np.array(dmp['alpha'], dtype=complex)
print("relaxation residual^2 at dumped alpha:", relax(al))
for a in [1.0, -1.0, 0.5, 2.0, 1j]:
    v = al.copy(); v[dmp['free_idx'][0]] = a
    print(f"  relaxation at a2 = {a}: {relax(v):.3e}")
# kernel-frame block of array 0
rows, cols, ent = closure(dmp, 0, isfree)
rmats = {}
for (l, k), r in rows.items():
    pass
# recover row matrices by recomputation
rowmat = {}
for (r2, c, i, jj) in ent:
    pass
print("closure array 0:", len(rows), "x", len(cols))
# terms of the last generator on array 0 at the kernel frame
for i, T in enumerate(gens):
    k0 = [h for (l, h) in T if l == 0 and isfree(h @ np.conj(K).T) is False and frame_of(0, h) == kf]
    if len(k0) >= 2:
        h0, h1 = k0[0], k0[1]
        u = K @ np.conj(h0).T @ h1 @ np.conj(K).T     # cell of col0 = K h0^-1 h1 = u K
        print(f"g{i}: {len(k0)} terms in the kernel frame; u = K h0^-1 h1 K^-1 free: {isfree(u)}")
        P = np.eye(u.shape[0], dtype=complex); o = 0
        while True:
            P = P @ u; o += 1
            if np.allclose(P, np.eye(u.shape[0])) or o > 64: break
        print(f"   order of u: {o}; scalar: {np.allclose(u, u[0,0]*np.eye(u.shape[0]))}; eigenvalues of u (f = 1): {np.round(sorted(np.linalg.eigvals(u), key=lambda z: (z.real, z.imag)), 3)}")
        # centrality in the free group: test against a few generators of B (Paulis X_e, Z_e, and an address map)
        n = u.shape[0]
        def X(z):
            M = np.zeros((n, n), dtype=complex)
            for x in range(n): M[x ^ z, x] = 1
            return M
        def Z(z):
            return np.diag([(-1) ** bin(x & z).count('1') for x in range(n)]).astype(complex)
        comm = all(np.allclose(u @ g, g @ u) for z in range(1, n) for g in (X(z), Z(z)))
        print(f"   u commutes with all Paulis: {comm}")
