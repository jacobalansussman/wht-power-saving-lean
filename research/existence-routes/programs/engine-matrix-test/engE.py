#!/usr/bin/env python3
"""engE.py (e3-engineE): the engine transfer form and the bordered rank identity (Theorem E10), checked on engine words.

Model (theory-guide hb.py, w-bounds section 1): block m, N = 2^m.  Letters are N x N matrices:
  dir z = ca I + cb X^z (paid; ca=(1+i)/2, cb=(1-i)/2), shifts X^z, address maps (GL(m,2) permutations),
  phase letters P_{z,c} = diag(i^{(z.x + c) mod 2}) (free).  K = dir e_1 ... dir e_m.
A word acts on registers; in a representation rho (here rho_f(g) = g^{(x)f}, the f-column action) the state is an
(R d) x (R d) matrix (register r = row block r).  The word is Lambda_0 D_1 Lambda_1 ... D_M Lambda_M with free products
Lambda_k and moves D_i = I + (rho(t_i) - 1) on register r_i.  With rho(t_i) - 1 = P_i Q_i^T (rank k_i):
  alpha_i = E_out^T Lambda_M..Lambda_i (e_{r_i} (x) P_i),  beta_i = (e_{r_i} (x) Q_i)^T Lambda_{i-1}..Lambda_0 E_in,
  G_{ji} = (e_{r_j} (x) Q_j)^T Lambda_{j-1}..Lambda_i (e_{r_i} (x) P_i) (i < j),  T_0 = E_out^T Lambda_M..Lambda_0 E_in
(moves deleted).  Transfer form: rho(total) = T_0 + A (I - G)^{-1} B.  Theorem E10: for every lambda,
  rank N(lambda) = sum_i k_i + rank(rho(total) - lambda),   N(lambda) = [[T_0 - lambda, A], [-B, I - G]].
"""
import itertools, sys
import numpy as np

class Eng:
    def __init__(s, m):
        s.m = m; N = s.N = 1 << m
        I = np.eye(N, dtype=complex); s.I = I
        def X(a):
            M = np.zeros((N, N), dtype=complex)
            for x in range(N): M[x ^ a, x] = 1
            return M
        def Zp(b): return np.diag([(-1) ** (bin(b & x).count("1") & 1) for x in range(N)]).astype(complex)
        s.X, s.Zp = X, Zp
        ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
        s.dir = {z: ca * I + cb * X(z) for z in range(1, N)}
        s.K = I.copy()
        for j in range(m): s.K = s.K @ s.dir[1 << j]
        s.free = {}
        for z in range(1, N): s.free["X%d" % z] = X(z)
        for z in range(N):
            for c in (0, 1):
                s.free["P%d_%d" % (z, c)] = np.diag([1j ** ((bin(z & x).count("1") + c) & 1) for x in range(N)])
        # a few address maps (random invertible) are added on demand
    def gl(s, cols):
        N = s.N; img = []
        for x in range(N):
            y = 0
            for k in range(s.m):
                if (x >> k) & 1: y ^= cols[k]
            img.append(y)
        assert len(set(img)) == N
        M = np.zeros((N, N), dtype=complex)
        for x in range(N): M[img[x], x] = 1
        return M
    def rand_gl(s, rng):
        while True:
            cols = [rng.randrange(1, s.N) for _ in range(s.m)]
            img = set()
            for x in range(s.N):
                y = 0
                for k in range(s.m):
                    if (x >> k) & 1: y ^= cols[k]
                img.add(y)
            if len(img) == s.N: return s.gl(cols)

def kpow(g, f):
    out = np.ones((1, 1), dtype=complex)
    for _ in range(f): out = np.kron(out, g)
    return out

def rank(M, tol=1e-8):
    if M.size == 0: return 0
    sv = np.linalg.svd(M, compute_uv=False)
    return int((sv > tol * max(1.0, sv[0])).sum())

def run_word(word, R, d, rep):
    """word: list of ('move', r, g) | ('free', r, g) | ('mix', r, [(s, mu), ...]) with g an N x N matrix.
    Returns the segments: list of Lambda matrices (Rd x Rd) and the list of moves (r, rho(g))."""
    Lams, moves = [], []
    S = np.eye(R * d, dtype=complex)
    def blk(r): return slice(r * d, (r + 1) * d)
    for st in word:
        if st[0] == "move":
            Lams.append(S); moves.append((st[1], rep(st[2]))); S = np.eye(R * d, dtype=complex)
        elif st[0] == "free":
            r = st[1]; S = S.copy(); S[blk(r)] = rep(st[2]) @ S[blk(r)]
        else:
            r, terms = st[1], st[2]
            new = sum(mu * S[blk(t)] for t, mu in terms) if terms else np.zeros((d, R * d), dtype=complex)
            S = S.copy(); S[blk(r)] = new
    Lams.append(S)
    return Lams, moves

def transfer(word, R, W, d, rep, outs, ins=None):
    if ins is None: ins = list(range(W))
    Lams, moves = run_word(word, R, d, rep)
    M = len(moves)
    def E(regs):
        Em = np.zeros((R * d, len(regs) * d), dtype=complex)
        for c, r in enumerate(regs): Em[r * d:(r + 1) * d, c * d:(c + 1) * d] = np.eye(d)
        return Em
    Ein, Eout = E(ins), E(outs)
    # total with moves
    tot = Lams[0]
    for i in range(M):
        r, g = moves[i]
        D = np.eye(R * d, dtype=complex); D[r * d:(r + 1) * d, r * d:(r + 1) * d] = g
        tot = Lams[i + 1] @ D @ tot
    total = Eout.T @ tot @ Ein
    # free product between positions
    def lprod(a, b):  # Lambda_b ... Lambda_a  (a <= b), product of free segments a..b
        P = np.eye(R * d, dtype=complex)
        for k in range(a, b + 1): P = Lams[k] @ P
        return P
    T0 = Eout.T @ lprod(0, M) @ Ein
    Ps, Qs, ks = [], [], []
    for r, g in moves:
        U, sv, Vh = np.linalg.svd(g - np.eye(d))
        k = int((sv > 1e-9).sum()); ks.append(k)
        e = np.zeros((R * d, d), dtype=complex); e[r * d:(r + 1) * d] = np.eye(d)
        Ps.append(e @ (U[:, :k] * sv[:k])); Qs.append(e @ Vh[:k].T)
    off = np.concatenate([[0], np.cumsum(ks)]).astype(int); Kt = int(off[-1])
    A = np.zeros((W * d, Kt), dtype=complex); B = np.zeros((Kt, W * d), dtype=complex); G = np.zeros((Kt, Kt), dtype=complex)
    for i in range(M):
        A[:, off[i]:off[i + 1]] = Eout.T @ lprod(i + 1, M) @ Ps[i]
        B[off[i]:off[i + 1], :] = Qs[i].T @ lprod(0, i) @ Ein
        for j in range(i + 1, M):
            G[off[j]:off[j + 1], off[i]:off[i + 1]] = Qs[j].T @ lprod(i + 1, j) @ Ps[i]
    # transfer form check
    tf = T0 + A @ np.linalg.solve(np.eye(Kt) - G, B) if Kt else T0
    return dict(total=total, T0=T0, A=A, B=B, G=G, ks=ks, tf_err=float(np.abs(tf - total).max()))

def bordered(tr, lam):
    T0, A, B, G = tr["T0"], tr["A"], tr["B"], tr["G"]
    Kt = G.shape[0]; n0 = T0.shape[0]
    N = np.zeros((n0 + Kt, n0 + Kt), dtype=complex)
    N[:n0, :n0] = T0 - lam * np.eye(n0); N[:n0, n0:] = A; N[n0:, :n0] = -B; N[n0:, n0:] = np.eye(Kt) - G
    return N
