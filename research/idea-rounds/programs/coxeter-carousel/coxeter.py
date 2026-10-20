"""Coxeter step checks (track wildcard, wht44/ideas, 16 Oct 2026).

Part A (frame game, exact over F_2): in Sp(2m,2) acting on Pauli labels (x|z), let h_j swap x_j<->z_j
(a paid move along direction j) and sigma cyclically relabel the m directions (a free address map).
The Coxeter step is c = sigma o h_1 (apply h_1, then sigma).  We check, for m = 2..8:
  * c^m = h_full (Hadamard on every direction), exactly;
  * c^(2m) = identity;
  * the frames F_j = (c^j)^-1 (L_Z) walk a geodesic L_Z -> L_X in m steps and come back in m more
    (a Petrie polygon of length 2m through the start and the target).
Part B (real arrays): with OpenAI's kernel C = (1/2)[[1+i,1-i],[1-i,1+i]] and with H, on arrays of
2^(m f) entries (m directions of f bits), (sigma o K_1)^m equals K on every direction, for m = 2,3 and
f = 1..3, to machine precision; (sigma o C_1)^(2m) is the free all-bits flip, (sigma o H_1)^(2m) = I.
"""
import itertools
import numpy as np

def sp_h(m, j):
    M = np.eye(2 * m, dtype=np.uint8)
    M[[j, m + j]] = M[[m + j, j]]
    return M

def sp_sigma(m):
    # Pauli X_j -> X_{j+1}, Z_j -> Z_{j+1}
    M = np.zeros((2 * m, 2 * m), dtype=np.uint8)
    for j in range(m):
        M[(j + 1) % m, j] = 1
        M[m + (j + 1) % m, m + j] = 1
    return M

def mul(A, B):
    return (A.astype(np.int64) @ B.astype(np.int64) % 2).astype(np.uint8)

def inv(A):
    n = A.shape[0]
    M = np.concatenate([A.copy(), np.eye(n, dtype=np.uint8)], axis=1)
    r = 0
    for col in range(n):
        piv = next(i for i in range(r, n) if M[i, col])
        M[[r, piv]] = M[[piv, r]]
        for i in range(n):
            if i != r and M[i, col]:
                M[i] ^= M[r]
        r += 1
    return M[:, n:]

def rank2(A):
    M = A.copy() % 2
    r = 0
    rows, cols = M.shape
    for col in range(cols):
        piv = [i for i in range(r, rows) if M[i, col]]
        if not piv:
            continue
        M[[r, piv[0]]] = M[[piv[0], r]]
        for i in range(rows):
            if i != r and M[i, col]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r

def lag_dist(L1, L2, m):
    # distance between Lagrangians (column spans) = m - dim(L1 cap L2) = rank([L1 L2]) - m
    return rank2(np.concatenate([L1, L2], axis=1)) - m

def part_a(out):
    for m in range(2, 9):
        h1 = sp_h(m, 0)
        hfull = np.eye(2 * m, dtype=np.uint8)
        for j in range(m):
            hfull = mul(sp_h(m, j), hfull)
        c = mul(sp_sigma(m), h1)
        P = np.eye(2 * m, dtype=np.uint8)
        powers = [P]
        for k in range(2 * m):
            P = mul(c, P)
            powers.append(P)
        LZ = np.zeros((2 * m, m), dtype=np.uint8)
        LZ[m:, :] = np.eye(m, dtype=np.uint8)
        LX = np.zeros((2 * m, m), dtype=np.uint8)
        LX[:m, :] = np.eye(m, dtype=np.uint8)
        dz, dx = [], []
        for k in range(2 * m + 1):
            F = mul(inv(powers[k]), LZ)
            dz.append(lag_dist(F, LZ, m))
            dx.append(lag_dist(F, LX, m))
        ok_m = np.array_equal(powers[m], hfull)
        ok_2m = np.array_equal(powers[2 * m], np.eye(2 * m, dtype=np.uint8))
        out.append(f"m={m}: c^m == h_full: {ok_m}; c^(2m) == I: {ok_2m}; "
                   f"d(F_j, L_Z) for j=0..2m: {dz}; d(F_j, L_X): {dx}")

def apply_dir(arr, K, j, m, f):
    """apply K^(x)f along direction j; arr has shape (2^f,)*m"""
    Kf = K
    for _ in range(f - 1):
        Kf = np.kron(Kf, K)
    return np.moveaxis(np.tensordot(Kf, arr, axes=([1], [j])), 0, j)

def sigma_arr(arr, m):
    # direction j -> direction j+1 (axis j of the input becomes axis j+1 of the output)
    return np.moveaxis(arr, list(range(m)), [(j + 1) % m for j in range(m)])

def part_b(out):
    H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    C = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    rng = np.random.default_rng(1)
    for name, K in (("H", H), ("C", C)):
        for m in (2, 3):
            for f in (1, 2, 3):
                shape = (2 ** f,) * m
                x = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
                y = x.copy()
                for _ in range(m):
                    y = sigma_arr(apply_dir(y, K, 0, m, f), m)
                t = x.copy()
                for j in range(m):
                    t = apply_dir(t, K, j, m, f)
                err_m = np.max(np.abs(y - t))
                z = y.copy()
                for _ in range(m):
                    z = sigma_arr(apply_dir(z, K, 0, m, f), m)
                ref = x.copy()
                if name == "C":
                    for j in range(m):
                        ref = apply_dir(ref, X, j, m, f)
                err_2m = np.max(np.abs(z - ref))
                out.append(f"{name}: m={m} f={f}: |(sigma K_1)^m x - K^(all) x| = {err_m:.1e}; "
                           f"|(sigma K_1)^(2m) x - {'Xflip' if name == 'C' else 'id'} x| = {err_2m:.1e}")

if __name__ == "__main__":
    out = ["# Part A: frame game Sp(2m,2), exact"]
    part_a(out)
    out.append("# Part B: real arrays, max abs error")
    part_b(out)
    txt = "\n".join(out)
    print(txt)
    open("out_coxeter.txt", "w").write(txt + "\n")
