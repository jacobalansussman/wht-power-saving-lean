"""openai_objects.py -- exact calibration of the objects in 08-reflection.tex, and the margin accounting.

Run: python3 openai_objects.py   (integer numpy, a few seconds)

Part A [calibration]: for d = 1..4 build C, Z, X, D_d, T_d, Gamma, J, M = Gamma(x)Gamma, L, S, A = S L and check
  every identity the paper states: D_d^2 = 0, Gamma^2 = I, nnz(Gamma) = (1 + d/2) N, J Gamma J^-1 = Gamma^T,
  A symmetric, nnz(L) = (d+2) w, (Gamma(x)I) L = I - M (so ker A = ker(I - M)), dim ker = w/2, the parity block form
  of M, and T_{m+1} = principal pivot transform of U (x) T_m on its bottom half (the exchange step of Lemma 8.1).
Part B [the converter]: with every 2x2 pair priced pi, the chain gives d w <= pi (d/2 + 4) w.  At pi = 2 (the LDU
  price) this is a tie up to 8w: every law except the dense pair is lossless at first order.  Threshold d* = 8 pi/(2-pi).
Part C [sparsity bound]: the feedback law gives p(G) <= nnz(G) for every invertible G (cell [[N, I],[I, 0]]).
  Check the cell identity numerically and report h / p(M) = nnz(Gamma) / (2 p(Gamma)) for OpenAI's Gamma.
Part D [ceilings]: zero-overhead ceilings as functions of pi: fixed-matrix (denominator ln w = (2d+2) ln 2) and the
  hypothetical uniform one (denominator ln(2d+2)).
"""
import numpy as np
from math import log

I2 = np.eye(2, dtype=np.int64)
C = np.array([[0, 0], [1, 0]], dtype=np.int64)
Z = np.array([[1, 0], [0, -1]], dtype=np.int64)
X = np.array([[0, 1], [1, 0]], dtype=np.int64)
U = np.array([[1, 0], [1, 1]], dtype=np.int64)


def kron(*ms):
    out = np.array([[1]], dtype=np.int64)
    for m in ms:
        out = np.kron(out, m)
    return out


def D(m):
    if m == 1:
        return C.copy()
    Dm = D(m - 1)
    n = Dm.shape[0]
    return np.block([[Dm, np.zeros((n, n), dtype=np.int64)], [np.eye(n, dtype=np.int64), -Dm]])


def rank_modp(Mat, p=2147483629):
    A = (Mat % p).astype(np.int64)
    r, c = A.shape
    rank = 0
    for col in range(c):
        piv = None
        for i in range(rank, r):
            if A[i, col] % p:
                piv = i
                break
        if piv is None:
            continue
        A[[rank, piv]] = A[[piv, rank]]
        inv = pow(int(A[rank, col]), p - 2, p)
        A[rank] = (A[rank] * inv) % p
        f = A[:, col].copy()
        f[rank] = 0
        nz = np.nonzero(f)[0]
        if len(nz):
            A[nz] = (A[nz] - (f[nz, None] * A[rank][None, :]) % p) % p
        rank += 1
        if rank == r:
            break
    return rank


def partA():
    print("Part A: OpenAI's reflection objects, exact integer checks")
    for d in range(1, 5):
        N = 2 ** (d + 1)
        w = N * N
        Dd = D(d)
        n = Dd.shape[0]
        In = np.eye(n, dtype=np.int64)
        Gam = np.block([[Dd, In], [In, -Dd]])
        ok = {}
        ok['D^2=0'] = not np.any(Dd @ Dd)
        ok['Gam^2=I'] = np.array_equal(Gam @ Gam, np.eye(N, dtype=np.int64))
        ok['nnz(Gam)=(1+d/2)N'] = (np.count_nonzero(Gam) * 2 == (2 + d) * N)
        J = kron(*[np.linalg.matrix_power(Z, i) @ X for i in range(d + 1)])
        ok['J Gam J^-1 = Gam^T'] = np.array_equal(J @ Gam @ J.T, Gam.T) and np.array_equal(J @ J.T, np.eye(N, dtype=np.int64))
        IN = np.eye(N, dtype=np.int64)
        M = np.kron(Gam, Gam)
        L = np.kron(Gam, IN) - np.kron(IN, Gam)
        S = np.kron(J, J)
        A = S @ L
        ok['A symmetric'] = np.array_equal(A, A.T)
        ok['A zero diagonal'] = not np.any(np.diag(A))
        ok['nnz(L)=(d+2)w'] = np.count_nonzero(L) == (d + 2) * w
        ok['(Gam x I) L = I - M'] = np.array_equal(np.kron(Gam, IN) @ L, np.eye(w, dtype=np.int64) - M)
        h = np.count_nonzero(A) // 2
        ok['h <= (d/2+1)w'] = 2 * h <= (d + 2) * w
        if w <= 256:
            rk = rank_modp(A)
            ok['rank A = w/2 (dim ker = w/2)'] = (rk == w // 2)
        # parity block form: M anticommutes with Pi (x) I, Pi = Z^(x)(d+1)
        Pi = kron(*([Z] * (d + 1)))
        ok['M anticommutes with Pi(x)I'] = np.array_equal(np.kron(Pi, IN) @ M, -M @ np.kron(Pi, IN))
        bad = [k for k, v in ok.items() if not v]
        print(f"  d={d}  N={N}  w={w}  h={h}  h/w={h/w:.3f}  (d/2+1)={d/2+1:.3f}  all {len(ok)} checks pass: {not bad}  {bad}")
    # T recursion = principal pivot transform of U (x) T_m on its bottom half
    for m in range(1, 5):
        Tm = np.eye(2 ** m, dtype=np.int64) + D(m)
        Tm_inv = np.eye(2 ** m, dtype=np.int64) - D(m)
        Tm1 = np.eye(2 ** (m + 1), dtype=np.int64) + D(m + 1)
        K0 = np.kron(U, Tm)               # [[Tm, 0], [Tm, Tm]]
        n = 2 ** m
        Ab, Bb, Cb, Db = K0[:n, :n], K0[n:, :n], K0[:n, n:], K0[n:, n:]
        Dinv = Tm_inv
        ppt = np.block([[Ab - Cb @ Dinv @ Bb, Cb @ Dinv], [-Dinv @ Bb, Dinv]])
        sgn = np.block([[np.eye(n, dtype=np.int64), np.zeros((n, n), dtype=np.int64)],
                        [np.zeros((n, n), dtype=np.int64), -np.eye(n, dtype=np.int64)]])
        same = np.array_equal(sgn @ ppt @ sgn, Tm1)
        nnzU = np.count_nonzero(np.kron(U, Tm)) - 2 * n
        print(f"  m={m}: T_(m+1) = diag(I,-I) ppt(U(x)T_m, bottom) diag(I,-I): {same};  off-diagonal nnz: "
              f"U^(x)(m+1) {3**(m+1)-2**(m+1)}, T_(m+1) {np.count_nonzero(Tm1) - 2**(m+1)}, price (m+1)2^m = {(m+1)*2**m}")


def partB():
    print("\nPart B: the converter.  All 2x2 pairs priced pi (dense pairs of K(eps), Cayley-cell scalar pairs,")
    print("        P_beta/Q_beta of the graph pivot, F and F^-1).  Chain: d w <= pi (d/2 + 4) w.")
    for pi, name in [(2.0, "LDU price 2"), (1.5, "dense-pair lemma 3/2"), (1.25, "5/4"), (1.0, "1"), (0.5, "floor 1/2 (U ~ B D B)")]:
        thr = 8 * pi / (2 - pi) if pi < 2 else float('inf')
        print(f"  pi = {pi:<5} ({name:<24}): margin per w = d(1 - pi/2) - 4 pi;  contradiction iff d > {thr:.2f}")
    print("  OpenAI's own bookkeeping (pairs at 3/2 only inside K(eps), 2 elsewhere): d w <= (3d/4 + 7.5) w, d > 30.")
    print("  With pi = 2 everywhere: d w <= (d + 8) w -- a tie up to 8w.  So the whole margin is (2 - pi) per pair of A.")


def partC():
    print("\nPart C: sparsity bound from the feedback law.  G invertible, D' generic diagonal, N = D'^-1 G - I invertible;")
    print("        cell K = [[N, I],[I, 0]] (= swap o rectangular shear, price <= nnz(N)) closes at s = 1 to I + N.")
    rng = np.random.default_rng(1)
    for trial in range(3):
        n = 6
        G = rng.integers(-2, 3, size=(n, n)) * (rng.random((n, n)) < 0.4) + np.diag(rng.integers(1, 3, size=n))
        Dp = np.diag(rng.integers(3, 7, size=n)).astype(float)
        Nn = np.linalg.inv(Dp) @ G - np.eye(n)
        K = np.block([[Nn, np.eye(n)], [np.eye(n), np.zeros((n, n))]])
        Mb, Cb, Bb, Db = K[:n, :n], K[:n, n:], K[n:, :n], K[n:, n:]
        H1 = Mb + Cb @ np.linalg.inv(np.eye(n) - Db) @ Bb
        print(f"  random G #{trial}: closure equals I + N: {np.allclose(H1, np.eye(n) + Nn)};  D' (I+N) == G: "
              f"{np.allclose(Dp @ (np.eye(n) + Nn), G)};  det N != 0: {abs(np.linalg.det(Nn)) > 1e-9};  nnz(N) = {np.count_nonzero(np.round(Nn, 12))}, nnz(G) = {np.count_nonzero(G)}")
    print("  Consequence: p(Gamma) <= nnz(Gamma); with A = S(Gamma(x)I - I(x)Gamma), h = N nnz(Gamma) >= N p(Gamma) = p(M)/2.")
    for d in [10, 30, 60, 100]:
        print(f"    OpenAI's Gamma at d={d}: p(Gamma)/(N) = {d/2}, nnz(Gamma)/N = {1+d/2}, h/p(M) = {(1+d/2)/d:.4f} (floor 1/2)")


def partD():
    print("\nPart D: zero-overhead ceilings, ratio r(d) = pi/2 + 4 pi/d")
    for pi in [1.5, 1.25, 1.0, 0.75, 0.5]:
        best_f = (0, None)
        best_u = (0, None)
        for d in range(2, 5000):
            r = pi / 2 + 4 * pi / d
            if r >= 1:
                continue
            sf = log(1 / r) / ((2 * d + 2) * log(2))
            su = log(1 / r) / log(2 * d + 2)
            if sf > best_f[0]:
                best_f = (sf, d)
            if su > best_u[0]:
                best_u = (su, d)
        print(f"  pi = {pi:<5}: fixed-matrix max s = {best_f[0]:.4e} at d = {best_f[1]};  uniform (hypothetical) max s = {best_u[0]:.4e} at d = {best_u[1]}")
    print("  (OpenAI's constants, pi = 3/2 inside K only: quantify/out_ceilings.txt gives 1.58e-3 fixed, 4.0e-2 uniform.)")


if __name__ == "__main__":
    partA()
    partB()
    partC()
    partD()
