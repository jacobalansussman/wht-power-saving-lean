"""check_model.py -- o3-close: cross-check the symplectic model (sp.py) against the 16x16 matrices of the exact engine
(o3-build/bm.py, fgeo.py letters).  For random words h = D b_k ... D b_1 b_0 (engine letters), compare
  depth(h) = log2 rowsupp(h)          with dist(sigma(h) L0, L0)
  dist_K(h) = log2 rowsupp(K h^-1)    with dist(sigma(h) Lambda_K, L0)
and sigma(letter) read off by conjugating Paulis."""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, HERE)
import numpy as np
from fgeo import Geo
import sp

def pauli(a, c):
    n = 16
    M = np.zeros((n, n), dtype=complex)
    for x in range(n):
        M[x ^ a, x] = (-1) ** bin(c & x).count('1')
    return M

def sigma_of(U):
    """symplectic matrix (tuple of images of basis vectors) of the Clifford U: w -> U P_w U^dag"""
    Ud = np.conj(U).T
    imgs = []
    for j in range(8):
        a, c = (1 << j, 0) if j < 4 else (0, 1 << (j - 4))
        P = U @ pauli(a, c) @ Ud
        col0 = P[:, 0]
        a2 = int(np.argmax(np.abs(col0)))
        ph = col0[a2]
        Z = pauli(a2, 0).T @ P / ph       # X_a2^-1 P / phase = Z_c2
        d = np.diag(Z)
        assert np.allclose(np.abs(Z - np.diag(d)), 0, atol=1e-9), "not a Pauli"
        c2 = 0
        for k in range(4):
            if abs(d[1 << k] - (-1)) < 1e-9:
                c2 |= 1 << k
        imgs.append(a2 | (c2 << 4))
    return tuple(imgs)

if __name__ == '__main__':
    geo = Geo(4, seed=7); G = geo.G
    rng = random.Random(7)
    assert sigma_of(G.D1) == sp.TAU, (sigma_of(G.D1), sp.TAU)
    assert sigma_of(G.K) == sp.KSP
    bad = 0; n = 0
    for trial in range(400):
        k = rng.randrange(1, 7)
        h = geo.rand_free()
        S = sigma_of(h)
        assert sp.img(S, sp.L0) == sp.L0
        for _ in range(k):
            b = geo.rand_free() if rng.random() < 0.8 else G.I
            h = G.mul(G.D1, b, h)
            S = sp.mat_mul(sp.TAU, sp.mat_mul(sigma_of(b), S))
        assert S == sigma_of(h)
        d_mat = geo.d0(h); k_mat = geo.dK(h)
        mu = sp.img(S, sp.L0); ka = sp.img(S, sp.LK)
        d_sp = sp.dist(mu, sp.L0); k_sp = sp.dist(ka, sp.L0)
        n += 1
        if (d_mat, k_mat) != (d_sp, k_sp):
            bad += 1
            print("MISMATCH", trial, (d_mat, k_mat), (d_sp, k_sp))
    print(f"checked {n} random words (1-6 moves): mismatches {bad}")
