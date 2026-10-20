# Independent check of (i) the f-column realisation of every binary quadratic form (point evaluation, shifts trivial,
# K -> (-1)^tr Q) and (ii) Lemma 1's transfer form on random programs (rectangular mixing, scratch, random move order).
import numpy as np, itertools
rng = np.random.default_rng(12345)
def nz(z, y):  # n(z) = #{c : z . y_c = 1}
    return sum(bin(z & yc).count("1") & 1 for yc in y)
def realise(q, m):
    y = []
    for (j, k) in q:
        if j == k: y += [1 << j, 1 << j]
        else: y += [1 << j, 1 << k, (1 << j) | (1 << k), 1 << j, 1 << j, 1 << k, 1 << k]
    return y
def Qval(q, z):
    v = 0
    for (j, k) in q: v ^= (z >> j & 1) & (z >> k & 1)
    return v
bad = 0; ntest = 0
for m in (3, 4, 5):
    pairs = [(j, k) for j in range(m) for k in range(j, m)]
    for _ in range(300):
        q = [p for p in pairs if rng.random() < 0.5]
        y = realise(q, m)
        tr = sum(1 for (j, k) in q if j == k) % 2
        for z in range(1, 1 << m):
            ntest += 1
            if nz(z, y) % 4 != 2 * Qval(q, z): bad += 1
            if (2 * nz(z, y)) % 4 != 0: bad += 1         # shift X^z = t_z^2 trivial
        Kval = 1j ** (sum(nz(1 << j, y) for j in range(m)) % 4)
        if abs(Kval - (-1) ** tr) > 1e-12: bad += 1
print("realisation: tests", ntest, "failures", bad)

# Lemma 1 on random programs: registers change size between moves.
maxerr = 0
for trial in range(200):
    m = 4; W = int(rng.integers(1, 4)); M = int(rng.integers(1, 9))
    sizes = [W] + [int(rng.integers(W, W + 4)) for _ in range(M)] + [W]
    Lam = [rng.standard_normal((sizes[i + 1], sizes[i])) + 1j * rng.standard_normal((sizes[i + 1], sizes[i])) for i in range(M + 1)]
    reg = [int(rng.integers(0, sizes[i + 1])) for i in range(M)]   # move i acts on register reg[i] of space i+1
    dirs = [int(rng.integers(1, 16)) for _ in range(M)]
    def prod(lo, hi):  # Lam[hi-1] ... Lam[lo] (identity parts of the moves in between)
        P = np.eye(sizes[lo])
        for t in range(lo, hi): P = Lam[t] @ P
        return P
    # move i (0-based) sits between Lam[i] and Lam[i+1], in space i+1
    al = np.array([prod(i + 1, M + 1)[:, reg[i]] for i in range(M)]).T     # W x M
    be = np.array([prod(0, i + 1)[reg[i], :] for i in range(M)])           # M x W
    G = np.zeros((M, M), complex)
    for i in range(M):
        for j in range(i + 1, M):
            G[j, i] = prod(i + 1, j + 1)[reg[j], reg[i]]
    Gam = 2 * G
    P0 = prod(0, M + 1)
    pairs = [(j, k) for j in range(m) for k in range(j, m)]
    for _ in range(20):
        q = [p for p in pairs if rng.random() < 0.5]
        cvals = [(-1) ** Qval(q, z) - 1 for z in dirs]
        P = np.eye(W)
        P = Lam[0] @ np.eye(W)
        for i in range(M):
            Dm = np.eye(sizes[i + 1], dtype=complex); Dm[reg[i], reg[i]] += cvals[i]
            P = Lam[i + 1] @ (Dm @ P)
        D = [i for i in range(M) if cvals[i] != 0]
        if D:
            F = P0 - 2 * al[:, D] @ np.linalg.solve(np.eye(len(D)) + Gam[np.ix_(D, D)], be[D, :])
        else:
            F = P0
        maxerr = max(maxerr, np.abs(F - P).max() / max(1, np.abs(P).max()))
print("Lemma 1 transfer form: 200 random programs x 20 forms, max relative error %.2e" % maxerr)
