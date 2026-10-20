#!/usr/bin/env python3
"""modgame.py: the one-sided (column) trivially pinned component, with the FULL free algebra, as a submodule game.
Frame module C^n (n = #frames), P = Stab(L0) acting by permutations; the space is (C^n)^W (one frame vector per input).
  V_0 = delta_{L0} (x) C^W;  move i: pick u_i in V_{i-1}; V_i = V_{i-1} + C[P] (M(t) - I) u_i.
The column condition M(Z) delta_{L0} = delta_{K L0} (x) I_W holds for some word with M moves iff some choice gives
delta_{K L0} (x) e_l in V_M for every l (outputs may be any element of V_M; links gamma in C[P] give any
u_i in V_{i-1}).  This program: generic choices, then a local search over the choices (coordinates of u_i in an
orthonormal basis of V_{i-1}) minimising the distance of the W targets to V_M.
usage: modgame.py m W M [restarts] [seed]"""
import sys, time
import numpy as np
from scipy.optimize import least_squares
from sp import Sp

def build(m):
    S = Sp(m); L0, F = S.frames(); idx = {L: i for i, L in enumerate(F)}; n = len(F)
    def perm_of(g): return np.array([idx[frozenset(S.apply(g, x) for x in L)] for L in F])
    gens = [perm_of(g) for g in S.P_gens()]
    t = perm_of(S.transvection(1))
    Kb = [((1 << j) if j < m else ((1 << j) ^ (1 << (j - m)))) for j in range(S.n)]
    return n, idx[L0], idx[frozenset(S.apply(Kb, x) for x in L0)], t, gens

def closure_inc(Q, gens, n, W, tol=1e-9):
    """NOT USED (tolerance drift: spurious dimensions). orthonormal basis ... (incremental:
    Gram-Schmidt with re-orthogonalisation; each new basis vector is pushed through every generator once)."""
    N = n * W
    B = np.zeros((N, 0)); queue = [Q[:, k] for k in range(Q.shape[1])]
    def add(v):
        nonlocal B
        for _ in range(2):
            v = v - B @ (B.T @ v)
        nv = np.linalg.norm(v)
        if nv > tol:
            v = v / nv; B = np.hstack([B, v[:, None]]); return v
        return None
    def act1(p, x):
        y = np.empty_like(x)
        for l in range(W): y[l * n:(l + 1) * n][p] = x[l * n:(l + 1) * n]
        return y
    scale = max([np.linalg.norm(q) for q in queue] + [1.0])
    queue = [q / scale for q in queue]
    while queue:
        v = add(queue.pop())
        if v is None: continue
        for p in gens: queue.append(act1(p, v))
    return B
def closure(Q, gens, n, W, tol=1e-9):
    def act(p, X):  # M(p) on each of the W blocks; X is (n W) x k
        Y = np.zeros_like(X)
        for l in range(W): Y[l * n:(l + 1) * n][p] = X[l * n:(l + 1) * n]
        return Y
    B, _ = np.linalg.qr(Q) if Q.shape[1] else (Q, None)
    B = B[:, np.linalg.svd(Q, compute_uv=False) > tol] if Q.shape[1] else B
    B = orth(Q, tol)
    while True:
        new = np.hstack([B] + [act(p, B) for p in gens])
        B2 = orth(new, tol)
        if B2.shape[1] == B.shape[1]: return B2
        B = B2

def orth(X, tol=1e-9):
    if X.shape[1] == 0: return X
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    return U[:, s > tol * max(1.0, s[0])]

def run(m, W, M, restarts=4, seed=0):
    n, i0, iK, t, gens = build(m)
    tinv = np.argsort(t)
    def Tm(x):  # (M(t) - I) on each block
        y = np.empty_like(x)
        for l in range(W): y[l * n:(l + 1) * n] = x[l * n:(l + 1) * n][tinv] - x[l * n:(l + 1) * n]
        return y
    V0 = np.zeros((n * W, W))
    for l in range(W): V0[l * n + i0, l] = 1
    T = np.zeros((n * W, W))
    for l in range(W): T[l * n + iK, l] = 1
    def play(coefs):
        V = V0; dims = [V.shape[1]]
        for i in range(M):
            c = coefs[i][:V.shape[1]] if len(coefs[i]) >= V.shape[1] else np.concatenate([coefs[i], np.zeros(V.shape[1] - len(coefs[i]))])
            u = V @ c
            w = Tm(u)
            V = closure(np.hstack([V, w[:, None]]), gens, n, W)
            dims.append(V.shape[1])
        R = T - V @ (V.T @ T)
        return np.linalg.norm(R), dims, V
    rng = np.random.default_rng(seed)
    dmax = n * W
    best = None
    for r in range(restarts):
        coefs = [rng.normal(size=dmax) for _ in range(M)]
        res, dims, V = play(coefs)
        if best is None or res < best[0]: best = (res, dims)
    return best, n

if __name__ == "__main__":
    m, W, M = map(int, sys.argv[1:4]); rs = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    t0 = time.time()
    (res, dims), n = run(m, W, M, rs)
    print("column-module game m=%d W=%d M=%d (frames %d, target needs dim >= %d): generic choices: dims %s, distance of targets to V_M = %.2e %s [%.1fs]"
          % (m, W, M, n, W * 2 ** (m * (m + 1) // 2), dims, res, "TARGETS REACHED" if res < 1e-8 else "", time.time() - t0), flush=True)
