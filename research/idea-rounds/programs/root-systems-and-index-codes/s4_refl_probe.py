"""s4_refl_probe.py -- numerical probe of the no-scratch orthogonal reflection form (evidence only, calibrated).

Form: W registers R^W, gates t = 0..M-1 in time order, gate t on direction z_t is the reflection
rho_t = I - 2 u_t u_t^T/|u_t|^2, active for the quadratic form Q iff Q(z_t) = 1. Requirement: for every Q,
rho_{M-1}^{q} ... rho_0^{q} = D_Q = diag((-1)^<Q,T_l>). (No scratch, invertible free maps; any such program is a
multi-target slice program, hence a no-phase unit by rm-theory Theorem B(a).)

Levenberg-Marquardt on all quadratic forms, analytic Jacobian via prefix/suffix products.
Usage: python3 s4_refl_probe.py MODE [seed]   MODE in {pos, neg, probe, sub7}
"""
import sys, itertools, time
import numpy as np

def setup(m):
    dirs = list(range(1, 2 ** m))
    zbits = np.array([[(z >> (m - 1 - i)) & 1 for i in range(m)] for z in dirs], dtype=np.int64)
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    Qs = np.array(list(itertools.product((0, 1), repeat=len(pairs))), dtype=np.int64)
    # Q(z) = sum_{i<=j} q_ij z_i z_j
    zz = np.array([[zb[i] * zb[j] for (i, j) in pairs] for zb in zbits], dtype=np.int64)  # 31 x npairs
    QZ = (Qs @ zz.T) % 2           # nQ x ndirs
    return dirs, zbits, pairs, Qs, QZ

def chars(Qs, pairs, T):
    v = np.array([T[i, j] % 2 for (i, j) in pairs], dtype=np.int64)
    return (Qs @ v) % 2

def products(U, act):
    """U: M x W gate vectors; act: nQ x M (0/1). Returns prefix list R[t] (product of gates < t) and P."""
    nQ, M = act.shape; W = U.shape[1]
    R = np.empty((M + 1, nQ, W, W)); R[0] = np.eye(W)
    for t in range(M):
        u = U[t] / np.linalg.norm(U[t]); rho = np.eye(W) - 2 * np.outer(u, u)
        A = act[:, t].astype(bool)
        R[t + 1] = R[t]
        R[t + 1][A] = np.einsum('ij,qjk->qik', rho, R[t][A])
    return R

def residual_jac(U, act, D):
    nQ, M = act.shape; W = U.shape[1]
    R = products(U, act)
    P = R[M]
    r = (P - D).reshape(nQ, W * W)
    # suffix S[t] = product of gates > t
    S = np.empty((M + 1, nQ, W, W)); S[M] = np.eye(W)
    rhos = []
    for t in range(M):
        u = U[t] / np.linalg.norm(U[t]); rhos.append(np.eye(W) - 2 * np.outer(u, u))
    for t in range(M - 1, -1, -1):
        A = act[:, t].astype(bool)
        S[t] = S[t + 1]
        S[t][A] = np.einsum('qij,jk->qik', S[t + 1][A], rhos[t])
    J = np.zeros((nQ * W * W, M * W))
    for t in range(M):
        A = act[:, t].astype(bool)
        if not A.any(): continue
        u = U[t]; n = u @ u
        for k in range(W):
            ek = np.zeros(W); ek[k] = 1
            drho = -2 * (np.outer(ek, u) + np.outer(u, ek)) / n + 4 * u[k] * np.outer(u, u) / n ** 2
            # S[t+1] drho R[t]
            G = np.einsum('qij,jk,qkl->qil', S[t + 1][A], drho, R[t][A])
            col = np.zeros((nQ, W, W)); col[A] = G
            J[:, t * W + k] = col.reshape(-1)
    return r.reshape(-1), J

def lm(U0, act, D, iters=400, tol=1e-24, verbose=False):
    U = U0.copy(); lam = 1e-2
    r, J = residual_jac(U, act, D); f = r @ r
    for it in range(iters):
        g = J.T @ r; Hm = J.T @ J
        while True:
            step = np.linalg.solve(Hm + lam * (np.diag(np.diag(Hm)) + 1e-12 * np.eye(len(g))), -g)
            Un = U + step.reshape(U.shape)
            rn, Jn = residual_jac(Un, act, D); fn = rn @ rn
            if fn < f:
                U, r, J, f = Un, rn, Jn, fn; lam = max(lam / 3, 1e-12); break
            lam *= 4
            if lam > 1e12: return U, f, it
        if verbose and it % 20 == 0: print("   it %d f %.6g lam %.1e" % (it, f, lam), flush=True)
        if f < tol: return U, f, it
    return U, f, iters

def rand_target(m, rng, need_nonalt=True):
    while True:
        Z = rng.integers(0, 2, size=(m, m))
        T = (Z @ Z.T) % 2
        # invertible over F_2 and congruent to I (nonzero diagonal)
        Mm = T.copy(); r = 0
        for c in range(m):
            p = next((i for i in range(r, m) if Mm[i, c]), None)
            if p is None: continue
            Mm[[r, p]] = Mm[[p, r]]
            for i in range(m):
                if i != r and Mm[i, c]: Mm[i] ^= Mm[r]
            r += 1
        if r == m and (not need_nonalt or T.diagonal().any()):
            return T, Z

def frame_dirs(Z, m):
    # columns of Z are the frame vectors z_i with T = sum z_i z_i^T
    return [int("".join(str(int(b)) for b in Z[:, i]), 2) for i in range(m)]

def run(mode, seed):
    rng = np.random.default_rng(seed)
    if mode in ("pos", "neg", "probe"):
        m = 4; dirs, zbits, pairs, Qs, QZ = setup(m); W = 4
        if mode == "pos":
            # 4 random frames (16 gates): the plain program exists
            Ts = []; gl = []
            for l in range(W):
                while True:
                    T, Z = rand_target(m, rng); fd = frame_dirs(Z, m)
                    if 0 not in fd and len(set(fd)) == m: break
                Ts.append(T); gl += fd
            order = rng.permutation(len(gl)); gl = [gl[i] for i in order]
        elif mode == "neg":
            # 13 random directions, targets with sum T_l = sum of the 13 z z^T (det condition), 4 in 13 (excluded by L24)
            gl = list(rng.choice(dirs, size=13, replace=True))
            S = sum(np.outer(zbits[z - 1], zbits[z - 1]) for z in gl) % 2
            while True:
                Ts = [rand_target(m, rng)[0] for _ in range(3)]
                T4 = (S + sum(Ts)) % 2
                if T4.diagonal().any() and round(abs(np.linalg.det(T4))) % 2 == 1 and \
                        len({tuple(t.flatten()) for t in Ts + [T4]}) == 4:
                    Ts.append(T4); break
        else:
            # 4 in 15, every direction once, distinct targets with sum 0 (det condition)
            while True:
                Ts = [rand_target(m, rng)[0] for _ in range(3)]
                T4 = sum(Ts) % 2
                if T4.diagonal().any() and round(abs(np.linalg.det(T4))) % 2 == 1 and \
                        len({tuple(t.flatten()) for t in Ts + [T4]}) == 4:
                    Ts.append(T4); break
            gl = list(rng.permutation(dirs))
        act = np.stack([QZ[:, z - 1] for z in gl], axis=1)
        D = np.zeros((len(Qs), W, W))
        for l, T in enumerate(Ts):
            D[:, l, l] = 1 - 2 * chars(Qs, pairs, T)
        best = None
        for start in range(4):
            U0 = rng.normal(size=(len(gl), W))
            U, f, it = lm(U0, act, D)
            print("%s seed %d start %d: M = %d gates, residual %.6g after %d iterations" % (mode, seed, start, len(gl), f, it), flush=True)
            best = f if best is None else min(best, f)
        print("%s seed %d BEST %.6g" % (mode, seed, best), flush=True)

if __name__ == "__main__":
    mode = sys.argv[1]; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    t0 = time.time(); run(mode, seed); print("time %.1f s" % (time.time() - t0))
