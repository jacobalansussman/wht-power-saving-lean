"""vmat.py -- v-o3-close: INDEPENDENT 16x16 engine matrices (no bm.py / fgeo.py) and an independent sigma map.

Letters: D = (1+i)/2 I + (1-i)/2 X_{e1}  (X_{e1} flips address bit 0);  K = D_{e1} D_{e2} D_{e3} D_{e4};
free letters = monomial Cliffords: affine address permutation y -> A y + c times a diagonal i^{q(y)},
q(y) = a.y + 2 sum_{j<k} B_jk y_j y_k + 2 z.y, times a scalar i^s.
depth(M) = log2 of the largest row support; dist_K(M) = depth(K M^-1).
sigma(U): conjugation of Paulis X^x Z^z, as a map on ints v = (x << 4) | z (vsp.py convention).
Run: python3 -I vmat.py NWORDS SEED   (checks the vsp model against the matrices)."""
import sys, os, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vsp

n = 16
I16 = np.eye(n, dtype=complex)

def Xmat(a):
    M = np.zeros((n, n), dtype=complex)
    for y in range(n):
        M[y ^ a, y] = 1
    return M

def Zmat(c):
    return np.diag([(-1) ** vsp.pc(c & y) for y in range(n)]).astype(complex)

def Dz(a):
    return (1 + 1j) / 2 * I16 + (1 - 1j) / 2 * Xmat(a)

D = Dz(1)
K = Dz(1) @ Dz(2) @ Dz(4) @ Dz(8)

def rand_gl(rng):
    while True:
        cols = [rng.randrange(1, 16) for _ in range(4)]
        if len(vsp.spanset(cols)) == 16:
            return cols

def rand_free(rng):
    cols = rand_gl(rng); c = rng.randrange(16)
    perm = []
    for y in range(n):
        v = c
        for j in range(4):
            if y >> j & 1:
                v ^= cols[j]
        perm.append(v)
    a = [rng.randrange(4) for _ in range(4)]
    B = {(j, k): rng.randrange(2) for j in range(4) for k in range(j + 1, 4)}
    z = rng.randrange(16); s = rng.randrange(4)
    M = np.zeros((n, n), dtype=complex)
    for y in range(n):
        bits = [(y >> j) & 1 for j in range(4)]
        q = sum(a[j] * bits[j] for j in range(4)) + 2 * sum(B[(j, k)] * bits[j] * bits[k] for (j, k) in B) \
            + 2 * vsp.pc(z & y) + s
        M[perm[y], y] = 1j ** (q % 4)
    return M

def rowsupp_max(M):
    return int((np.abs(M) > 1e-9).sum(axis=1).max())

def depth(M):
    r = rowsupp_max(M)
    d = int(round(np.log2(r)))
    assert 1 << d == r, r
    return d

def distK(M):
    return depth(K @ np.conj(M).T)

def is_monomial(M):
    nz = np.abs(M) > 1e-9
    return bool((nz.sum(axis=1) == 1).all() and (nz.sum(axis=0) == 1).all())

def sigma(U):
    """returns a function-table: list img[v] for v in 0..255 (linear map on F_2^8)"""
    Ud = np.conj(U).T
    basis_img = []
    for j in range(8):
        v = 1 << j
        P = U @ Xmat(vsp.X(v)) @ Zmat(vsp.Z(v)) @ Ud
        col0 = P[:, 0]
        x2 = int(np.argmax(np.abs(col0)))
        R = Xmat(x2) @ P                         # should be phase * Z^{z2}
        dg = np.diag(R)
        assert np.allclose(R, np.diag(dg), atol=1e-9)
        ph = dg[0]
        z2 = 0
        for k in range(4):
            val = dg[1 << k] / ph
            if abs(val + 1) < 1e-9:
                z2 |= 1 << k
            else:
                assert abs(val - 1) < 1e-9
        basis_img.append((x2 << 4) | z2)
    tab = []
    for v in range(256):
        w = 0
        for j in range(8):
            if v >> j & 1:
                w ^= basis_img[j]
        tab.append(w)
    return tab

def sg_apply(tab, S):
    return frozenset(tab[v] for v in S)

if __name__ == '__main__':
    NW = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 11)
    sD = sigma(D); sK = sigma(K)
    assert all(sD[v] == vsp.tau(v) for v in range(256)), "sigma(D) != tau"
    assert all(sK[v] == vsp.sigK(v) for v in range(256)), "sigma(K) != model"
    # free letters keep L0
    for _ in range(50):
        b = rand_free(rng); assert sg_apply(sigma(b), vsp.L0) == vsp.L0
    bad = 0; hist = {}
    for w in range(NW):
        k = rng.randrange(1, 8)
        h = rand_free(rng)
        for _ in range(k):
            b = rand_free(rng) if rng.random() < 0.85 else I16
            h = D @ b @ h
        S = sigma(h)
        mu = sg_apply(S, vsp.L0); ka = sg_apply(S, vsp.LK)
        pm = (depth(h), distK(h)); pmod = (vsp.dist(mu, vsp.L0), vsp.dist(ka, vsp.L0))
        hist[pm] = hist.get(pm, 0) + 1
        if pm != pmod:
            bad += 1; print("MISMATCH", w, pm, pmod)
    print(f"words {NW} (1-7 moves, own matrices, own sigma): model mismatches {bad}")
    print("(depth, dist_K) histogram:", sorted(hist.items()))
