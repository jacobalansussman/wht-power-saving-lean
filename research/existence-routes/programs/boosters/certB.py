# e3-boost certB.py: profiles of the label-permutation certificates A and B (w-bounds 5.4) at blocks m = 2..5.
# A: V = Q^(nonzero Pauli labels), Y_U = span of the unit vectors of U's nonzero labels.  delta_A(d) = 2^m - 2^(m-d).
# B: A modulo N = span{s_U - s_V : frames U, V} (s_U = sum of U's unit vectors).  w-bounds 5.4 found 3, 5, 6 at m = 3.
# Conjecture checked here: delta_B(d) = 2^m - 2^(m-d) - 1 at every m.  Ranks mod two primes (w.h.p. = rank over Q;
# at m <= 4 N is spanned by ALL frames, at m = 5 by random frames until the rank stops growing).
import sys, random
import numpy as np

def omega(a, b, m):
    M = (1 << m) - 1
    ax, az, bx, bz = a >> m, a & M, b >> m, b & M
    return bin((ax & bz) ^ (az & bx)).count("1") & 1

def frame_span(gens):
    S = {0}
    for g in gens:
        S |= {s ^ g for s in S}
    return frozenset(S)

def L_d(m, d):  # x-part e_1..e_d, z-part e_{d+1}..e_m: a frame at distance d from L0 = {(0|z)}
    return frame_span([(1 << i) << m for i in range(d)] + [1 << j for j in range(d, m)])

def all_frames(m):
    L0 = L_d(m, 0); seen = {L0}; todo = [L0]; n = 1 << (2 * m)
    while todo:
        U = todo.pop();
        for v in range(1, n):
            W = frozenset(u ^ v if omega(u, v, m) else u for u in U)
            if W not in seen: seen.add(W); todo.append(W)
    return list(seen)

def rand_frame(m, steps=40):
    U = L_d(m, 0); n = 1 << (2 * m)
    for _ in range(steps):
        v = random.randrange(1, n); U = frozenset(u ^ v if omega(u, v, m) else u for u in U)
    return U

class RankMod:
    def __init__(self, n, p): self.n, self.p, self.B, self.piv = n, p, np.zeros((0, n), np.int64), []
    def add(self, R):
        p = self.p; R = np.array(R, np.int64) % p
        if self.piv:
            R = (R - (R[:, self.piv] @ self.B) % p) % p   # B is reduced: identity on pivot columns
        for row in R:
            for (c, b) in zip(self.piv, self.B):
                if row[c]: row = (row - row[c] * b) % p
            nz = np.nonzero(row)[0]
            if len(nz) == 0: continue
            c = nz[0]; row = (row * pow(int(row[c]), p - 2, p)) % p
            if len(self.B):
                col = self.B[:, c].copy()
                self.B = (self.B - np.outer(col, row)) % p
            self.B = np.vstack([self.B, row]); self.piv.append(c)
        return len(self.piv)
    def copy(self):
        o = RankMod(self.n, self.p); o.B = self.B.copy(); o.piv = list(self.piv); return o

def vec(U, n):
    v = np.zeros(n - 1, np.int64)
    for u in U:
        if u: v[u - 1] = 1
    return v

def run(m, p, nrand=None):
    n = 1 << (2 * m)
    frames = all_frames(m) if nrand is None else [rand_frame(m) for _ in range(nrand)]
    s0 = vec(L_d(m, 0), n)
    N = RankMod(n - 1, p)
    rows = [vec(U, n) - s0 for U in frames]
    rk = 0
    for i in range(0, len(rows), 400):
        rk = N.add(rows[i:i + 400])
    Y0 = [np.eye(n - 1, dtype=np.int64)[u - 1] for u in L_d(m, 0) if u]
    NY = N.copy(); base = NY.add(Y0)
    out = []
    for d in range(1, m + 1):
        T = NY.copy(); Yd = [np.eye(n - 1, dtype=np.int64)[u - 1] for u in L_d(m, d) if u]
        out.append(T.add(Yd) - base)
    return len(frames), rk, base, out

if __name__ == "__main__":
    random.seed(1)
    for m in [2, 3, 4] + ([5] if len(sys.argv) > 1 else []):
        for p in (1000003, 999983):   # products < 2^40, sums of <= 1023 terms < 2^50: no int64 overflow
            nf, rkN, base, dB = run(m, p, None if m <= 4 else int(sys.argv[1]))
            dA = [2 ** m - 2 ** (m - d) for d in range(1, m + 1)]
            print(f"m={m} p={p}: frames used {nf}, dim N {rkN}, dim(N+Y0) {base}; delta_A {dA}; delta_B {dB}; "
                  f"delta_B == delta_A - 1: {dB == [a - 1 for a in dA]}", flush=True)
