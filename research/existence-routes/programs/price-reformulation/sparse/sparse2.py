#!/usr/bin/env python3
"""Cleaner version of sparse_annih.py: every real symmetric A with A E = 0 is A = Q Y Q^T (Q orthonormal basis of
E^perp, Y symmetric).  Reweighted L1 on the off-diagonal entries of A, normalization <Y, Y0> = 1 with Y0 random
positive definite (or Y0 = I), then: support -> generic element on that support -> exact rank / kernel check.
Usage: python3 sparse2.py d iters seed
"""
import sys, time
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import null_space
from scipy.sparse import lil_matrix, csr_matrix, hstack, vstack, identity
d = int(sys.argv[1]); iters = int(sys.argv[2]); seed = int(sys.argv[3])
C = np.array([[0., 0.], [1., 0.]]); Z = np.diag([1., -1.]); X = np.array([[0., 1.], [1., 0.]])
def kron(*ms):
    out = np.array([[1.]])
    for m in ms: out = np.kron(out, m)
    return out
D = C.copy()
for m in range(1, d): D = kron(C, np.eye(2 ** m)) + kron(Z, D)
Gam = kron(X, np.eye(2 ** d)) + kron(Z, D)
N = 2 ** (d + 1); w = N * N; r = w // 2
J = np.array([[1.]])
for i in range(1, d + 2): J = np.kron(J, np.linalg.matrix_power(Z, i - 1) @ X)
L = np.kron(Gam, np.eye(N)) - np.kron(np.eye(N), Gam)
A0 = np.kron(J, J) @ L
K = null_space(L); Q = null_space(K.T)       # Q: orthonormal basis of E^perp (w x r)
def pairs(A, tol):
    B = (np.abs(A) > tol) | (np.abs(A.T) > tol)
    return int(np.triu(B, 1).sum())
h0 = pairs(A0, 1e-9)
yidx = [(a, b) for a in range(r) for b in range(a, r)]; ny = len(yidx)
pidx = [(i, j) for i in range(w) for j in range(i + 1, w)]; npairs = len(pidx)
# linear map Y -> A_ij = sum_{a<=b} Y_ab * (Q_ia Q_jb + Q_ib Q_ja) (a<b) or Q_ia Q_ja (a=b)
M = np.zeros((npairs, ny))
for k, (a, b) in enumerate(yidx):
    if a == b:
        col = np.outer(Q[:, a], Q[:, a])
    else:
        col = np.outer(Q[:, a], Q[:, b]) + np.outer(Q[:, b], Q[:, a])
    M[:, k] = col[np.triu_indices(w, 1)]
g = np.random.default_rng(seed)
G = g.standard_normal((r, r)); Y0 = G @ G.T / r + np.eye(r)
nrow = np.array([Y0[a, b] * (1 if a == b else 2) for (a, b) in yidx])
print(f"d={d} w={w}: OpenAI h0={h0} ({h0/w:.3f} w); LP vars {ny}+{npairs}")
wts = np.ones(npairs); best = None
for it in range(iters):
    t0 = time.time()
    # variables [y (free), t (>=0)]; constraints  M y - t <= 0, -M y - t <= 0 ; nrow.y = 1
    Aub = np.block([[M, -np.eye(npairs)], [-M, -np.eye(npairs)]])
    bub = np.zeros(2 * npairs)
    Aeq = np.concatenate([nrow, np.zeros(npairs)])[None, :]
    c = np.concatenate([np.zeros(ny), wts])
    bounds = [(None, None)] * ny + [(0, None)] * npairs
    res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1.0], bounds=bounds, method='highs')
    if res.status != 0:
        print("  LP status", res.status, res.message); break
    y = res.x[:ny]; a = M @ y
    tol = 1e-7 * np.abs(a).max()
    supp = np.abs(a) > tol
    h = int(supp.sum())
    # generic element of {Y : A_ij(Y) = 0 off the support}
    Voff = M[~supp]
    V = null_space(Voff) if Voff.shape[0] else np.eye(ny)
    yy = V @ g.standard_normal(V.shape[1])
    Ymat = np.zeros((r, r))
    for k, (p_, q_) in enumerate(yidx): Ymat[p_, q_] = Ymat[q_, p_] = yy[k]
    rkY = np.linalg.matrix_rank(Ymat, tol=1e-9 * np.abs(Ymat).max())
    B = Q @ Ymat @ Q.T
    hb = pairs(B, 1e-9 * np.abs(B).max())
    valid = rkY == r
    print(f"  iter {it}: L1 support h={h} ({h/w:.3f} w); generic element: h={hb}, rank={rkY}/{r} valid={valid} dim(V)={V.shape[1]} [{time.time()-t0:.1f}s]")
    if valid and (best is None or hb < best): best = hb
    wts = 1.0 / (np.abs(a) / np.abs(a).max() + 1e-3)
print(f"RESULT d={d} seed={seed}: OpenAI h0={h0} ({h0/w:.3f} w); best valid h={best} ({(best/w) if best else float('nan'):.3f} w)")
