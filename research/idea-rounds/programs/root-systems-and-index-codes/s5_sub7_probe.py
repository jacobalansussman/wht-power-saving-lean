"""s5_sub7_probe.py -- the 7-array subspace problem (rmfam N14) in the no-scratch orthogonal reflection form.

Block 5, GF(32) natural form (as s3_subspace.py); arrays indexed by v in V \ 0, target T_v = mult by beta_v^2.
Registers: W = 7 (arrays only) or W = 8 (plus one clean register owed the identity).
LM on a random subset of the 32,768 quadratic forms (analytic Jacobian), then a check on all of them.
Usage: python3 s5_sub7_probe.py MODE W seed nstarts   MODE in {pos (plain 7 in 35), probe (7 in 31), ord (7 in 31,
       time order = the plain frames' directions grouped by array)}
Evidence only: a stall decides nothing; a hit must be checked on all forms and then exactly.
"""
import sys, itertools, time
import numpy as np

M_ = 5; Q_ = 32; N_ = 31; POLY = 0b100101
exp = [0] * (2 * N_); log = [None] * Q_
x = 1
for k in range(N_):
    exp[k] = x; log[x] = k; x <<= 1
    if x & Q_: x ^= POLY
for k in range(N_, 2 * N_): exp[k] = exp[k - N_]
def mul(a, b): return 0 if a == 0 or b == 0 else exp[log[a] + log[b]]
def tr(a):
    s = 0; y = a
    for _ in range(M_): s ^= y; y = mul(y, y)
    return s
def self_dual():
    for cand in itertools.combinations(range(1, Q_), M_):
        if all(tr(mul(a, b)) == (1 if a == b else 0) for a in cand for b in cand): return list(cand)
F = self_dual()
def coords(z): return np.array([tr(mul(z, f)) for f in F], dtype=np.int64)
pairs = [(i, j) for i in range(5) for j in range(i, 5)]
Qs = np.array(list(itertools.product((0, 1), repeat=15)), dtype=np.int64)
zz = np.array([[coords(z)[i] * coords(z)[j] for (i, j) in pairs] for z in range(1, 32)], dtype=np.int64)
QZ = (Qs @ zz.T) % 2   # 32768 x 31, column z-1
S = [0, 1, 2, 5, 11, 18, 19]
betas = [exp[k] for k in S]
def mmat(g): return np.array([[tr(mul(g, mul(fi, fj))) for fj in F] for fi in F], dtype=np.int64)
Ts = [mmat(mul(b, b)) for b in betas]
chi = np.stack([(Qs @ np.array([T[i, j] for (i, j) in pairs])) % 2 for T in Ts], axis=1)   # 32768 x 7

def products(U, act):
    nQ, M = act.shape; W = U.shape[1]
    R = np.empty((M + 1, nQ, W, W)); R[0] = np.eye(W)
    rhos = []
    for t in range(M):
        u = U[t] / np.linalg.norm(U[t]); rho = np.eye(W) - 2 * np.outer(u, u); rhos.append(rho)
        A = act[:, t].astype(bool)
        R[t + 1] = R[t]; R[t + 1][A] = np.matmul(rho, R[t][A])
    return R, rhos

def residual_jac(U, act, D, want_jac=True):
    nQ, M = act.shape; W = U.shape[1]
    R, rhos = products(U, act)
    r = (R[M] - D).reshape(-1)
    if not want_jac: return r, None
    Sx = np.empty((M + 1, nQ, W, W)); Sx[M] = np.eye(W)
    for t in range(M - 1, -1, -1):
        A = act[:, t].astype(bool)
        Sx[t] = Sx[t + 1]; Sx[t][A] = np.matmul(Sx[t + 1][A], rhos[t])
    J = np.zeros((nQ * W * W, M * W))
    for t in range(M):
        A = act[:, t].astype(bool)
        if not A.any(): continue
        u = U[t]; n = u @ u
        SA = Sx[t + 1][A]; RA = R[t][A]
        for k in range(W):
            ek = np.zeros(W); ek[k] = 1
            drho = -2 * (np.outer(ek, u) + np.outer(u, ek)) / n + 4 * u[k] * np.outer(u, u) / n ** 2
            G = np.matmul(np.matmul(SA, drho), RA)
            col = np.zeros((nQ, W, W)); col[A] = G
            J[:, t * W + k] = col.reshape(-1)
    return r, J

def lm(U0, act, D, iters=int(__import__("os").environ.get("ITERS", "150")), tol=1e-22):
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
        if f < tol: return U, f, it
    return U, f, iters

def run(mode, W, seed, nstarts, nsub=1024):
    rng = np.random.default_rng(seed)
    if mode == "pos":
        gl = [mul(b, f) for b in betas for f in F]          # the 7 plain frames, 35 gates
        gl = [gl[i] for i in rng.permutation(len(gl))]
    elif mode == "ord":
        # 7 in 31 with the time order of the directions grouped as in the plain frames (each direction once)
        seen = []
        for b in betas:
            for f in F:
                z = mul(b, f)
                if z not in seen: seen.append(z)
        rest = [z for z in range(1, 32) if z not in seen]
        gl = seen + rest
    else:
        gl = list(rng.permutation(np.arange(1, 32)))
    act_all = np.stack([QZ[:, z - 1] for z in gl], axis=1)
    D_all = np.zeros((len(Qs), W, W))
    for l in range(7): D_all[:, l, l] = 1 - 2 * chi[:, l]
    if W == 8: D_all[:, 7, 7] = 1
    sub = rng.choice(len(Qs), size=nsub, replace=False)
    act, D = act_all[sub], D_all[sub]
    for start in range(nstarts):
        t0 = time.time()
        U0 = rng.normal(size=(len(gl), W))
        U, f, it = lm(U0, act, D)
        rfull, _ = residual_jac(U, act_all, D_all, want_jac=False)
        print("%s W=%d seed %d start %d: M=%d, subset residual %.6g (%d its), ALL-forms residual %.6g, %.0f s"
              % (mode, W, seed, start, len(gl), f, it, rfull @ rfull, time.time() - t0), flush=True)

if __name__ == "__main__":
    mode = sys.argv[1]; W = int(sys.argv[2]); seed = int(sys.argv[3]); nstarts = int(sys.argv[4])
    run(mode, W, seed, nstarts)
