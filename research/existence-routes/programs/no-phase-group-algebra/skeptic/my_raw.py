#!/usr/bin/env python3
"""my_raw.py (v-e3-rmbuild): own raw-array checks (no hb.py, no e3-rmbuild code).
Block m, f digit columns: an array is a complex tensor of shape (2^m,)*f; a matrix letter acts on every factor.
  X(z)[x, x^z] = 1, dir(z) = (1+i)/2 I + (1-i)/2 X(z), K = dir(e_1)...dir(e_m), P_g[g x, x] = 1,
  c_z = (1-i)/2 I + (1+i)/2 X(z)  (a mixing of the register with its shifted copy: free).
Checks: (1) Lemma R: c_z dir(z) = Hadamard-basis sign (-1)^binom(s,2), s = #{c : z.y_c = 1}; (2) Fano: product of
the corrected moves over all nonzero z of F_2^3 is 1; (3) P_h d'_z P_h^-1 = d'_{hz}; (4) lifted words of
C[Sym_3 x| H] (U_3 = 2 in 7, a two-channel 2 in 6 with an address map) deliver K y_l exactly; uncorrected word fails.
Usage: my_raw.py fmax_m2 fmax_m3"""
import sys
import numpy as np

rng = np.random.default_rng(20261015)


def mats(m):
    N = 1 << m
    def X(z):
        A = np.zeros((N, N), complex)
        for x in range(N):
            A[x, x ^ z] = 1
        return A
    I = np.eye(N)
    D = {z: (1 + 1j) / 2 * I + (1 - 1j) / 2 * X(z) for z in range(N)}
    C = {z: (1 - 1j) / 2 * I + (1 + 1j) / 2 * X(z) for z in range(N)}
    K = I.astype(complex)
    for j in range(m):
        K = K @ D[1 << j]
    return X, D, C, K


def lin(m, cols):
    N = 1 << m; img = []
    for x in range(N):
        y = 0
        for k in range(m):
            if x >> k & 1:
                y ^= cols[k]
        img.append(y)
    P = np.zeros((N, N), complex)
    for x in range(N):
        P[img[x], x] = 1
    return P, img


def apc(z, v, X):
    """free element c_z = (1-i)/2 [1] + (1+i)/2 [X^z]: one scalar mixing of v with its shifted copy."""
    return (1 - 1j) / 2 * v + (1 + 1j) / 2 * ap(X(z), v)


def ap(M, v):
    f = v.ndim
    for ax in range(f):
        v = np.moveaxis(np.tensordot(M, v, axes=([1], [ax])), 0, ax)
    return v


def rand(m, f):
    sh = (1 << m,) * f
    return rng.standard_normal(sh) + 1j * rng.standard_normal(sh)


def hadamard(m):
    N = 1 << m
    return np.array([[(-1) ** bin(x & y).count("1") for y in range(N)] for x in range(N)], float)


def check_R(m, f):
    X, D, C, K = mats(m); N = 1 << m; Hd = hadamard(m)
    worst = 0
    for z in range(1, N):
        v = rand(m, f)
        lhs = apc(z, ap(D[z], v), X)
        # expected: in Hadamard basis multiply by (-1)^binom(s,2)
        w = ap(Hd, v)
        idx = np.indices((N,) * f)
        par = np.array([bin(y & z).count("1") & 1 for y in range(N)])
        s = sum(par[idx[c]] for c in range(f))
        sign = np.where((s * (s - 1) // 2) % 2 == 0, 1.0, -1.0)
        exp = ap(Hd / N, sign * w)
        worst = max(worst, np.max(np.abs(lhs - exp)) / np.max(np.abs(v)))
        # involution
        v2 = apc(z, ap(D[z], lhs), X)
        worst = max(worst, np.max(np.abs(v2 - v)) / np.max(np.abs(v)))
    return worst


def dprime(m, z, v):
    X, D, C, K = mats(m)
    return apc(z, ap(D[z], v), X)


def run_sym_word(m, word, W, R, f, corrected=True):
    """word over C[Sym x| H]: ('smove', r, z) | ('P', r, Pmatrix) | ('mix', r, [(s, coef), ...]).
    Lift: smove -> dir(z) then c_z (if corrected).  End: delivering registers multiplied by prod_j c_{e_j}^-1."""
    X, D, C, K = mats(m)
    sh = (1 << m,) * f
    regs = [np.zeros(sh, complex) for _ in range(R)]
    ys = [rand(m, f) for _ in range(W)]
    for l in range(W):
        regs[l] = ys[l].copy()
    for st in word:
        if st[0] == "smove":
            _, r, z = st
            regs[r] = ap(D[z], regs[r])
            if corrected:
                regs[r] = apc(z, regs[r], X)
        elif st[0] == "P":
            regs[st[1]] = ap(st[2], regs[st[1]])
        elif st[0] == "mix":
            new = np.zeros(sh, complex)
            for s_, c in st[2]:
                new = new + c * regs[s_]
            regs[st[1]] = new
    return regs, ys, (C, X, K)


def check_word(m, word, W, R, deliver, f, corrected=True):
    regs, ys, (C, X, K) = run_sym_word(m, word, W, R, f, corrected)
    err = 0
    for l, r in enumerate(deliver):
        out = regs[r]
        if corrected:
            for j in range(m):
                out = ap(X(1 << j), apc(1 << j, out, X))        # c^-1 = c X
        err = max(err, np.max(np.abs(out - ap(K, ys[l]))) / np.max(np.abs(ys[l])))
    return err


if __name__ == "__main__":
    f2, f3 = int(sys.argv[1]), int(sys.argv[2])
    f3min = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    for f in range(1, f2 + 1):
        print("Lemma R m=2 f=%d: max rel err %.2e" % (f, check_R(2, f))); sys.stdout.flush()
    for f in range(f3min, f3 + 1):
        print("Lemma R m=3 f=%d: max rel err %.2e" % (f, check_R(3, f))); sys.stdout.flush()
    m = 3
    X, D, C, K = mats(m)
    for f in range(f3min, f3 + 1):
        v = rand(m, f); w = v.copy()
        for z in range(1, 8):
            w = dprime(m, z, w)
        e1 = np.max(np.abs(w - v)) / np.max(np.abs(v))
        P, img = lin(m, [3, 6, 7]); Pi = np.linalg.inv(P)
        e2 = 0
        for z in range(1, 8):
            v = rand(m, f)
            a = ap(P, dprime(m, z, ap(Pi, v)))
            b = dprime(m, img[z], v)
            e2 = max(e2, np.max(np.abs(a - b)) / np.max(np.abs(v)))
        print("m=3 f=%d: Fano product of 7 corrected moves = 1: err %.2e; P_h d'_z P_h^-1 = d'_hz: err %.2e" % (f, e1, e2))
        sys.stdout.flush()
    # lifted words
    U3 = [("smove", 0, z) for z in (1, 2, 4)] + [("smove", 1, z) for z in (3, 5, 6, 7)]
    Pc, _ = lin(m, [2, 4, 1]); Pci = np.linalg.inv(Pc)
    T7 = [("P", 1, Pc), ("mix", 2, [(0, 1), (1, 1)]), ("mix", 3, [(0, 1), (1, -1)])]
    T7 += [("smove", r, z) for r in (2, 3) for z in (1, 2, 4)]
    T7 += [("mix", 0, [(2, .5), (3, .5)]), ("mix", 1, [(2, .5), (3, -.5)]), ("P", 1, Pci)]
    # two channels with a NON-orthogonal map (should fail on array 2)
    Pn, _ = lin(m, [3, 6, 7]); Pni = np.linalg.inv(Pn)
    T8 = [("P", 1, Pn)] + T7[1:-1] + [("P", 1, Pni)]
    for f in range(f3min, f3 + 1):
        print("f=%d: U_3 lifted err %.2e | U_3 uncorrected err %.2e | T7 lifted err %.2e | T8 (must fail) err %.2e"
              % (f, check_word(m, U3, 2, 2, [0, 1], f), check_word(m, U3, 2, 2, [0, 1], f, corrected=False),
                 check_word(m, T7, 2, 4, [0, 1], f), check_word(m, T8, 2, 4, [0, 1], f)))
        sys.stdout.flush()
