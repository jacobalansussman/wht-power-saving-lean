#!/usr/bin/env python3
"""Frobenius (angle-weighted) ceiling for frame families on the dual polar graph DSp(2m,2).

For a distance-transitive graph, every G-invariant PSD kernel K on vertices averages to
K_d = sum_i c_i omega_i(d) with c_i >= 0, omega_i the cosine sequence (spherical function) of eigenvalue theta_i.
phi_d := K_0 - K_d = sum_{i>=1} c_i (1 - omega_i(d)).  Hence phi_j/phi_1 <= max_{i>=1} (1-omega_i(j))/(1-omega_i(1)).

Intersection array of DSp(2m,2) (classical parameters, q = 2, e = 1):
  b_j = 2^(j+1) (2^(m-j) - 1),  c_j = 2^j - 1,  k = b_0,  a_j = k - b_j - c_j.
Eigenvalues: theta_i = 2(2^(m-i) - 1) - (2^i - 1)   (Theorem S' of c2-finite-group; checked here by the closing equation).
Exact rational arithmetic throughout.
"""
from fractions import Fraction as F
import sys

def arrays(m):
    b = [2 ** (j + 1) * (2 ** (m - j) - 1) for j in range(m + 1)]
    c = [2 ** j - 1 for j in range(m + 1)]
    k = b[0]
    a = [k - b[j] - c[j] for j in range(m + 1)]
    return k, a, b, c

def cosines(m, theta):
    k, a, b, c = arrays(m)
    w = [F(1), F(theta, k)]
    for j in range(1, m):
        # c_j w_{j-1} + a_j w_j + b_j w_{j+1} = theta w_j
        w.append((theta * w[j] - a[j] * w[j] - c[j] * w[j - 1]) / b[j])
    # closing equation at j = m (b_m = 0) must hold iff theta is an eigenvalue
    close = c[m] * w[m - 1] + a[m] * w[m] - theta * w[m]
    return w, close

def main(mmax=30):
    print("m  T_m=max_i (1-w_i(m))/(1-w_i(1))   argmax i   2-2^(1-m)   equal?   per-rank ceilings T_m(j), j=1..m (first 6)")
    for m in range(1, mmax + 1):
        thetas = [2 * (2 ** (m - i) - 1) - (2 ** i - 1) for i in range(m + 1)]
        W = []
        for i, th in enumerate(thetas):
            w, close = cosines(m, th)
            assert close == 0, (m, i, close)
            W.append(w)
        # multiplicity check not needed for the bound; only cosines enter
        Tj = []
        for j in range(1, m + 1):
            best, arg = None, None
            for i in range(1, m + 1):
                r = (1 - W[i][j]) / (1 - W[i][1])
                if best is None or r > best:
                    best, arg = r, i
            Tj.append((best, arg))
        T, arg = Tj[-1]
        target = 2 - F(2, 2 ** m)
        print(f"{m:2d}  {str(T):>24s} = {float(T):.6f}   i={arg}   {float(target):.6f}   {T == target}   "
              + ", ".join(f"{float(t):.4f}(i={a})" for t, a in Tj[:6]))
        if m <= 6:
            for i in range(m + 1):
                print("      i=%d theta=%d  omega_i(0..m) = %s" % (i, thetas[i], [str(x) for x in W[i]]))

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 30)

def closed_form_checks(mmax=40):
    """(a) omega_1(j) = (2^(m+1-j) + 2^(1-j) - 4)/(2^(m+1) - 2)  [proved in PROOF.md for all m];
       (b) omega_i(m) = (-1)^i 2^(-i(m+1-i)) for 0 <= i <= m (single-term opposition eigenvalue);
       (c) hence T_m = 2 - 2^(1-m), attained at i = 1 only;
       (d) multiplicities m_i = N / sum_j k_j omega_i(j)^2 are positive integers summing to N = prod (2^j + 1)."""
    for m in range(1, mmax + 1):
        k, a, b, c = arrays(m)
        thetas = [2 * (2 ** (m - i) - 1) - (2 ** i - 1) for i in range(m + 1)]
        kj = [1]
        for j in range(m): kj.append(kj[-1] * b[j] // c[j + 1])
        N = 1
        for j in range(1, m + 1): N *= 2 ** j + 1
        assert sum(kj) == N, (m, sum(kj), N)
        W = [cosines(m, th)[0] for th in thetas]
        for j in range(m + 1):
            assert W[1][j] == F(2 ** (m + 1 - j) + F(2, 2 ** j) - 4, 2 ** (m + 1) - 2), (m, j)
        for i in range(m + 1):
            assert W[i][m] == F((-1) ** i, 2 ** (i * (m + 1 - i))), (m, i, W[i][m])
        mults = []
        for i in range(m + 1):
            s = sum(kj[j] * W[i][j] ** 2 for j in range(m + 1))
            mi = F(N) / s
            assert mi.denominator == 1 and mi > 0, (m, i, mi)
            mults.append(int(mi))
        assert sum(mults) == N
        ratios = [(1 - W[i][m]) / (1 - W[i][1]) for i in range(1, m + 1)]
        T = max(ratios)
        assert T == 2 - F(2, 2 ** m) and ratios.index(T) == 0 and ratios.count(T) == 1, (m, ratios)
        if m <= 8 or m % 8 == 0:
            print(f"m={m}: N={N} sphere sizes {kj[:5]}{'...' if m > 4 else ''} multiplicities {mults if m <= 6 else mults[:4]} "
                  f"T_m={T} ratios i=1..: {[round(float(x), 5) for x in ratios[:5]]}")
    print(f"closed forms (a)-(d) verified exactly for m = 1..{mmax}")

if __name__ == "__main__":
    closed_form_checks(40)
