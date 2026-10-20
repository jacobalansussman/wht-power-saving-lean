#!/usr/bin/env python3
"""v_e3.py (skeptic v-e3-engineE): exact checks of Proposition E3 and Lemma E0 with Gaussian-integer matrices
(numpy complex with small integer entries: exact).  Own code.
 E3: 2^m K == Had P_K Had  (Had = +-1 Hadamard, P_K = diag(i^|x|)), and 2 S_j dir(e_j) S_j == (1+i) Had_j, so
     h = prod_j S_j dir(e_j) S_j (a product of letters) satisfies h P_K h^-1 = K.  Then, in the f-column action
     g -> g^(x)f, rank(rho(P_K) - lam) == rank(rho(K) - lam) for f = 1..3 (consequence, checked numerically).
 E0: dir z == g dir(e_1) g^-1 for an address map g with g e_1 = z, for every z != 0.
 E10 spot check (own word): a random engine word at m = 2, f = 1, 2, transfer form and bordered rank identity.
usage: v_e3.py"""
import itertools
import numpy as np

def setup(m):
    N = 1 << m
    def X(z):
        A = np.zeros((N, N), complex)
        for x in range(N): A[x ^ z, x] = 1
        return A
    def dir2(z): return (1 + 1j) * np.eye(N) + (1 - 1j) * X(z)     # 2 * dir z
    def Sj(j): return np.diag([1j ** ((x >> j) & 1) for x in range(N)])
    def Hadj(j):
        A = np.zeros((N, N), complex)
        for x in range(N):
            A[x, x] += (-1) ** ((x >> j) & 1); A[x ^ (1 << j), x] += 1
        return A
    PK = np.diag([1j ** bin(x).count("1") for x in range(N)])
    Had = np.array([[(-1) ** bin(x & y).count("1") for y in range(N)] for x in range(N)], complex)
    return N, X, dir2, Sj, Hadj, PK, Had

for m in (1, 2, 3, 4):
    N, X, dir2, Sj, Hadj, PK, Had = setup(m)
    K2 = np.eye(N, dtype=complex)
    for j in range(m): K2 = K2 @ dir2(1 << j)                       # 2^m K
    ok1 = np.array_equal(Had @ PK @ Had, K2)
    ok2 = all(np.array_equal(2 * Sj(j) @ (dir2(1 << j) / 2) @ Sj(j), (1 + 1j) * Hadj(j)) for j in range(m))
    # Hadj product = Had
    Hp = np.eye(N, dtype=complex)
    for j in range(m): Hp = Hp @ Hadj(j)
    ok3 = np.array_equal(Hp, Had)
    # Lemma E0
    def addr(A):
        def Ax(x):
            o = 0
            for j in range(m):
                if (x >> j) & 1: o ^= A[j]
            return o
        G = np.zeros((N, N), complex)
        for x in range(N): G[Ax(x), x] = 1
        return G
    okE0 = True
    for z in range(1, N):
        # an invertible A with A e_1 = z: complete z to a basis
        cols = [z]
        for e in range(m):
            cand = 1 << e
            span = {0}
            for c in cols: span |= {s ^ c for s in span}
            if cand not in span and len(cols) < m: cols.append(cand)
        G = addr(cols)
        okE0 &= np.array_equal(G @ dir2(1) @ np.linalg.inv(G), dir2(z))
    print("m=%d: Had P_K Had == 2^m K: %s; 2 S_j dir(e_j) S_j == (1+i) Had_j for all j: %s; prod Had_j == Had: %s; "
          "Lemma E0 (dir z = g dir e1 g^-1, all z != 0): %s" % (m, ok1, ok2, ok3, okE0))

# rank consequence in f-column actions
for m, fs in ((2, (1, 2, 3)), (3, (1, 2))):
    N, X, dir2, Sj, Hadj, PK, Had = setup(m)
    K = np.eye(N, dtype=complex)
    for j in range(m): K = K @ (dir2(1 << j) / 2)
    out = []
    for f in fs:
        rK, rP = K, PK
        for _ in range(f - 1): rK, rP = np.kron(rK, K), np.kron(rP, PK)
        same = all(np.linalg.matrix_rank(rK - lam * np.eye(len(rK)), tol=1e-8) == np.linalg.matrix_rank(rP - lam * np.eye(len(rP)), tol=1e-8)
                   for lam in (1, 1j, -1, -1j, 0.3 + 0.7j))
        out.append((f, same))
    print("m=%d: rank(rho_f(P_K) - lam) == rank(rho_f(K) - lam), lam in {1,i,-1,-i,0.3+0.7i}: %s" % (m, out))

# Theorem E10 spot check on an own random word: registers r = 0..R-1; inputs in registers 0..W-1; outputs W..2W-1? keep
# simple: W = 1, R = 3, M moves on random registers, free steps = random C[B] mixing matrices (combinations of free
# letters with random complex coefficients), evaluated in rho_f.
rng = np.random.default_rng(99)
m = 2; N, X, dir2, Sj, Hadj, PK, Had = setup(m)
def phase(z, c): return np.diag([1j ** ((bin(z & x).count("1") + c) % 2) for x in range(N)])
freeL = [np.eye(N, dtype=complex)] + [X(z) for z in range(1, N)] + [phase(z, c) for z in range(N) for c in (0, 1)]
P01 = np.zeros((N, N), complex)
for x in range(N): P01[((x & 1) << 1) | ((x >> 1) & 1), x] = 1   # swap address bits (address map)
freeL.append(P01)
t = dir2(1) / 2
for f in (1, 2):
    def rho(g):
        o = g
        for _ in range(f - 1): o = np.kron(o, g)
        return o
    d = N ** f; R = 3; M = 4; W = 1
    # Lambda_k: R x R blocks, each a random combination of 2 free letters
    def randLam():
        L = np.zeros((R * d, R * d), complex)
        for r in range(R):
            for s in range(R):
                for _ in range(2):
                    L[r * d:(r + 1) * d, s * d:(s + 1) * d] += (rng.normal() + 1j * rng.normal()) * rho(freeL[rng.integers(len(freeL))])
        return L
    Lam = [randLam() for _ in range(M + 1)]
    regs = rng.integers(R, size=M)
    def Dmove(r):
        D = np.eye(R * d, dtype=complex); D[r * d:(r + 1) * d, r * d:(r + 1) * d] = rho(t); return D
    Tot = Lam[0]
    for i in range(M): Tot = Lam[i + 1] @ Dmove(regs[i]) @ Tot
    Ein = np.zeros((R * d, d), complex); Ein[:d] = np.eye(d)
    Eout = Ein
    Z = Eout.T @ Tot @ Ein
    # transfer form pieces
    def prod(a, b):  # Lambda_a ... Lambda_b (a >= b)
        P = np.eye(R * d, dtype=complex)
        for k in range(b, a + 1): P = Lam[k] @ P
        return P
    def er(r):
        E = np.zeros((R * d, d), complex); E[r * d:(r + 1) * d] = np.eye(d); return E
    T0 = Eout.T @ prod(M, 0) @ Ein
    U, s, Vh = np.linalg.svd(rho(t) - np.eye(d)); k = int((s > 1e-10).sum())
    Pf = U[:, :k] * s[:k]; Qt = Vh[:k]                                   # rho(t) - 1 = Pf Qt
    alpha = [Eout.T @ prod(M, i + 1) @ er(regs[i]) for i in range(M)]
    beta = [er(regs[i]).T @ prod(i, 0) @ Ein for i in range(M)]
    At = np.hstack([alpha[i] @ Pf for i in range(M)])
    Bt = np.vstack([Qt @ beta[i] for i in range(M)])
    Gm = np.zeros((M * k, M * k), complex)
    for i in range(M):
        for j in range(i + 1, M):
            gam = er(regs[j]).T @ prod(j, i + 1) @ er(regs[i])
            Gm[j * k:(j + 1) * k, i * k:(i + 1) * k] = Qt @ gam @ Pf
    Zt = T0 + At @ np.linalg.solve(np.eye(M * k) - Gm, Bt)
    errs = np.abs(Zt - Z).max()
    rk = []
    for lam in (1, 1j, -1, 0.37 + 0.81j, 0):
        Nl = np.block([[T0 - lam * np.eye(d), At], [-Bt, np.eye(M * k) - Gm]])
        rk.append(np.linalg.matrix_rank(Nl, tol=1e-8) == M * k + np.linalg.matrix_rank(Z - lam * np.eye(d), tol=1e-8))
    print("E10 own random word m=2 f=%d (M=4 moves, 3 registers, k=%d): transfer form error %.1e; rank identity at 5 lambdas: %s"
          % (f, k, errs, rk))
