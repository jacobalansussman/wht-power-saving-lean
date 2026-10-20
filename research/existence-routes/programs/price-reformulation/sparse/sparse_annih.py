#!/usr/bin/env python3
"""Is OpenAI's 'cheap side' improvable?  Their upper bound p(R) <= (3/2) h + 2w counts the supported unordered
off-diagonal pairs h of a REAL SYMMETRIC matrix A with ker A = E (E = positive eigenspace of M = Gamma (x) Gamma).
They use A0 = (J (x) J)(Gamma (x) I - I (x) Gamma), h0 <= (d/2 + 1) w.  Here: search for a sparser symmetric A with
ker A = E exactly (reweighted L1 for the support, then a random element of the subspace on that support, then an
exact rank check).  d = 1, 2 (w = 16, 64), d = 3 (w = 256) if time allows.
Usage: python3 sparse_annih.py d [iters]
"""
import sys, time
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import null_space

d = int(sys.argv[1]); iters = int(sys.argv[2]) if len(sys.argv) > 2 else 6
mode = sys.argv[3] if len(sys.argv) > 3 else 'A0'; seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
C = np.array([[0., 0.], [1., 0.]]); Z = np.diag([1., -1.]); X = np.array([[0., 1.], [1., 0.]])
def kron(*ms):
    out = np.array([[1.]])
    for m in ms: out = np.kron(out, m)
    return out
D = C.copy()
for m in range(1, d):
    D = kron(C, np.eye(2 ** m)) + kron(Z, D)
Gam = kron(X, np.eye(2 ** d)) + kron(Z, D)
N = 2 ** (d + 1); w = N * N
assert np.allclose(Gam @ Gam, np.eye(N))
J = np.array([[1.]])
for i in range(1, d + 2):
    J = np.kron(J, np.linalg.matrix_power(Z, i - 1) @ X)
assert np.allclose(J @ Gam @ np.linalg.inv(J), Gam.T)
L = np.kron(Gam, np.eye(N)) - np.kron(np.eye(N), Gam)
S = np.kron(J, J)
A0 = S @ L
assert np.allclose(A0, A0.T)
K = null_space(L)            # basis of E, w x w/2
assert K.shape[1] == w // 2, K.shape
def pairs(A, tol=1e-9):
    B = (np.abs(A) > tol) | (np.abs(A.T) > tol)
    return int(np.triu(B, 1).sum()), int(np.diag(np.abs(A) > tol).sum())
h0, dg0 = pairs(A0)
print(f"d={d} N={N} w={w}: OpenAI A0 has h0={h0} off-diagonal pairs (h0/w={h0/w:.3f}; bound d/2+1={d/2+1}), diag nnz={dg0}, rank={np.linalg.matrix_rank(A0)}")
# variables: upper triangle incl. diagonal
idx = [(i, j) for i in range(w) for j in range(i, w)]
nv = len(idx)
pos = {p: k for k, p in enumerate(idx)}
# constraint A K = 0 : for each row r and each column c of K: sum_j A[r,j] K[j,c] = 0
rows = []
Aeq = np.zeros((w * K.shape[1] + 1, nv))
r = 0
for i in range(w):
    for c in range(K.shape[1]):
        for j in range(w):
            k = pos[(min(i, j), max(i, j))]
            Aeq[r, k] += K[j, c]
        r += 1
# normalization: mode A0: <A, A0>_F = <A0, A0>_F ; mode rand: <A, R>_F = 1 for R = P A_rand P, a random valid
# symmetric matrix (P = projector onto E^perp), so the LP is not steered towards A0
if mode == 'A0':
    T = A0; val = float((A0 * A0).sum())
else:
    g = np.random.default_rng(seed)
    P = np.eye(w) - K @ K.T
    G = g.standard_normal((w, w)); G = G + G.T
    T = P @ G @ P; val = 1.0
for (i, j), k in pos.items():
    Aeq[r, k] = T[i, j] * (1 if i == j else 2)
beq = np.zeros(Aeq.shape[0]); beq[-1] = val
print(f"  normalization mode={mode} seed={seed}")
# drop dependent equality rows (numerically) to help the LP
Q, R_, P = __import__('scipy.linalg', fromlist=['qr']).qr(Aeq.T, mode='economic', pivoting=True)
rk = int((np.abs(np.diag(R_)) > 1e-9 * abs(R_[0, 0])).sum())
keep = np.sort(P[:rk])
if (Aeq.shape[0] - 1) not in keep: keep = np.sort(np.append(keep[:-1], Aeq.shape[0] - 1))
Aeq_r, beq_r = Aeq[keep], beq[keep]
print(f"  LP: {nv} variables, {Aeq.shape[0]} equations, rank {rk}")
offd = np.array([1.0 if i != j else 0.0 for (i, j) in idx])
wts = offd.copy()
best = None
for it in range(iters):
    t0 = time.time()
    # x = xp - xn, minimize sum wts*(xp+xn)
    c = np.concatenate([wts, wts]) + 1e-12
    res = linprog(c, A_eq=np.hstack([Aeq_r, -Aeq_r]), b_eq=beq_r, bounds=(0, None), method='highs')
    if res.status != 0:
        print("  LP failed", res.message); break
    x = res.x[:nv] - res.x[nv:]
    A = np.zeros((w, w))
    for (i, j), k in pos.items():
        A[i, j] = A[j, i] = x[k]
    h, dg = pairs(A, 1e-7)
    # support -> subspace -> random element -> exact kernel check
    supp = [k for k, (i, j) in enumerate(idx) if abs(x[k]) > 1e-7 or i == j]
    sub = Aeq[:-1][:, supp]
    V = null_space(sub)
    rng = np.random.default_rng(it)
    y = V @ rng.standard_normal(V.shape[1])
    B = np.zeros((w, w))
    for kk, k in enumerate(supp):
        i, j = idx[k]; B[i, j] = B[j, i] = y[kk]
    hb, _ = pairs(B, 1e-9)
    rkB = np.linalg.matrix_rank(B, tol=1e-8 * np.abs(B).max())
    okker = np.allclose(B @ K, 0, atol=1e-8 * np.abs(B).max())
    valid = (rkB == w // 2) and okker
    print(f"  iter {it}: L1 solution has h={h} (h/w={h/w:.3f}); generic element on its support: h={hb}, rank={rkB} (need {w//2}), BK=0: {okker} -> valid={valid}  [{time.time()-t0:.1f}s]")
    if valid and (best is None or hb < best): best = hb
    wts = offd / (np.abs(x) + 1e-3)
print(f"RESULT d={d}: OpenAI h0={h0} ({h0/w:.3f} w); best valid sparse h found = {best} ({(best/w if best else float('nan')):.3f} w); ratio {((best/h0) if best else float('nan')):.3f}")
