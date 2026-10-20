#!/usr/bin/env python3
"""t_lemmaR.py (e3-rmbuild): exact check of Lemma R on raw arrays (small, Mac-safe).

Claim: c_z := ((1-i)[1] + (1+i)[X^z]) / 2  (a free combination of shifts) gives the corrected move
    d'_z := c_z * dir(z)   (formal product: ((1-i)[dir z] + (1+i)[X^z dir z]) / 2),
which at f digit columns acts in the X-basis as the sign (-1)^{binom(s,2)}, s = #{columns c : z.y_c = 1}.
Hence d'_z is an involution, the d'_z commute, P_h d'_z P_h^-1 = d'_{hz}, and at every character in every shift-sector
d'_z acts as (-1)^{Q(z)} (Q the sector's binary quadratic form): all sectors are the same problem.
Checked: m = 2: f = 1..5; m = 3: f = 1..3; every z.
Matrices are built at f columns as Kronecker powers of the f = 1 matrices; entries are Gaussian dyadic rationals,
compared with tolerance 1e-12 (they are exactly representable in binary floating point at these sizes).
"""
import itertools
import numpy as np


def mats(m):
    N = 1 << m
    I = np.eye(N, dtype=complex)
    def X(a):
        M = np.zeros((N, N), dtype=complex)
        for x in range(N):
            M[x ^ a, x] = 1
        return M
    d = {z: ((1 + 1j) * I + (1 - 1j) * X(z)) / 2 for z in range(1, N)}
    return I, X, d


def kronpow(M, f):
    R = np.array([[1.0 + 0j]])
    for _ in range(f):
        R = np.kron(R, M)
    return R


def main():
    bad = 0; checks = 0
    for m, fmax in ((2, 5), (3, 3)):
        N = 1 << m
        I, X, d = mats(m)
        H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        Hm = kronpow(H, m)
        for f in range(1, fmax + 1):
            Hf = kronpow(Hm, f)
            for z in range(1, N):
                dz = kronpow(d[z], f); xz = kronpow(X(z), f)
                dp = ((1 - 1j) * dz + (1 + 1j) * (xz @ dz)) / 2
                # X-basis
                D = Hf @ dp @ Hf.conj().T
                # expected sign: columns y_c in F_2^m, s = #{c : z.y_c = 1}
                exp = np.zeros(N ** f)
                for idx, ys in enumerate(itertools.product(range(N), repeat=f)):
                    s = sum(bin(z & y).count("1") & 1 for y in ys)
                    exp[idx] = (-1) ** (s * (s - 1) // 2)
                # Kronecker index order: first factor most significant -> itertools.product order matches
                ok = np.allclose(D, np.diag(exp), atol=1e-12)
                inv = np.allclose(dp @ dp, np.eye(N ** f), atol=1e-12)
                checks += 1
                if not (ok and inv):
                    bad += 1
                    print("FAIL", m, f, z, ok, inv)
        print("m=%d: checked f=1..%d, all z" % (m, fmax))
    print("checks", checks, "failures", bad)


if __name__ == "__main__":
    main()
