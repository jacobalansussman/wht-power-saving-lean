#!/usr/bin/env python3
"""sep_algebra.py (mw2-singer4, step 1): separation matrices over C[<eta>], eta of order 7 (a Singer cycle of F_2^3).

Elements of Z[i][Z7] are length-7 arrays of Gaussian integers (a + b i stored as two int arrays).  For a W x W matrix T
with such entries:
  * complex invertibility: det T(zeta^k) != 0 for the 7 characters k = 0..6 (k = 0 is the trivial representation);
  * 2-adic reduction: Z[i] -> F_2 (pi = 1 - i, a + bi -> a + b mod 2), then F_2[Z7] = F_2 x F_8 x F_8' by the CRT
    (t^7 - 1 = (t + 1)(t^3 + t + 1)(t^3 + t^2 + 1)); "singular mod 2" = det vanishes in at least one factor.
    The factor pattern is a subset of {F2, F8, F8'}: F2 = the trivial character mod 2 (augmentation parity), F8 / F8' =
    evaluation at a root alpha of t^3 + t + 1 / at alpha^-1 (root of t^3 + t^2 + 1).
  * exact inverse over Q(i)[Z7] (block-circulant solve with Fractions), its denominators.
Symmetries used for classes: multiplying rows/columns by units +-eta^k, +-i eta^k, row/column permutations, Galois
eta -> eta^a (a in (Z/7)^*; a = -1 swaps F8 and F8'), complex conjugation.
Matching with RM(1,3): a 7-bit word f (coefficients of det mod 2, cyclic order) extended by a parity bit lies in an RM(1,3)
coset determined by the syndrome (f(alpha), f(1)) in F8 x F2 [the cyclic [7,4] Hamming code is the ideal (t^3+t+1)];
so "det kills F8" <-> det mod 2 is a Hamming codeword, "det kills F2" <-> even weight.
usage: python3 sep_algebra.py  (writes out_sep_algebra.txt next to it)"""
import itertools, sys, os, cmath
from fractions import Fraction as Fr
from collections import Counter, defaultdict
import numpy as np

N = 7
ZETA = [cmath.exp(2j * cmath.pi * k / N) for k in range(N)]

# ---- F_8 = F_2[a]/(a^3 + a + 1), elements as 3-bit ints
def f8mul(x, y):
    r = 0
    for i in range(3):
        if (y >> i) & 1:
            r ^= x << i
    for d in (4, 3):
        if (r >> d) & 1:
            r ^= 0b1011 << (d - 3)
    return r
ALPHA = 0b010
APOW = [1]
for _ in range(N - 1):
    APOW.append(f8mul(APOW[-1], ALPHA))
assert f8mul(APOW[-1], ALPHA) == 1
def ev8(bitsf, power):
    """evaluate f in F_2[t]/(t^7-1) at alpha^power"""
    r = 0
    for j in range(N):
        if bitsf[j]:
            r ^= APOW[(j * power) % N]
    return r

# ---- group ring Z[i][Z7]: element = complex numpy array of length 7 (integer real/imag parts)
def mul(a, b):
    out = np.zeros(N, dtype=complex)
    for i in range(N):
        if a[i] != 0:
            out += a[i] * np.roll(b, i)
    return out

def det(T):
    W = len(T)
    if W == 1:
        return T[0][0]
    out = np.zeros(N, dtype=complex)
    for j in range(W):
        minor = [row[:j] + row[j + 1:] for row in T[1:]]
        t = mul(T[0][j], det(minor))
        out = out + t if j % 2 == 0 else out - t
    return out

def chars(a):
    return [sum(a[j] * ZETA[(j * k) % N] for j in range(N)) for k in range(N)]

def mod2(a):
    return [int(round(a[j].real) + round(a[j].imag)) % 2 for j in range(N)]

def residue_pattern(a):
    f = mod2(a)
    kill = []
    if sum(f) % 2 == 0:
        kill.append("F2")
    if ev8(f, 1) == 0:
        kill.append("F8")
    if ev8(f, 6) == 0:          # alpha^-1 = alpha^6
        kill.append("F8'")
    return tuple(kill), f

def complex_invertible(a, tol=1e-9):
    return all(abs(c) > tol for c in chars(a))

# ---- exact inverse over Q(i)[Z7] of a W x W matrix: solve the (2*7*W) real block system with Fractions
def to_real_block(T):
    W = len(T)
    n = N * W
    M = [[Fr(0)] * (2 * n) for _ in range(2 * n)]
    for r in range(W):
        for c in range(W):
            e = T[r][c]
            for j in range(N):
                if e[j] == 0:
                    continue
                re, im = int(round(e[j].real)), int(round(e[j].imag))
                for s in range(N):          # e * eta^s = sum_j e_j eta^{j+s}: column (c, s) -> row (r, j+s)
                    R0 = r * N + (j + s) % N; C0 = c * N + s
                    M[R0][C0] += re; M[R0][C0 + n] -= im
                    M[R0 + n][C0] += im; M[R0 + n][C0 + n] += re
    return M

def inverse_denominators(T):
    """exact inverse; returns (lcm of denominators, max 2-adic exponent of denominators) or None if singular over Q"""
    M = to_real_block(T)
    n2 = len(M)
    A = [row[:] + [Fr(int(i == j)) for j in range(n2)] for i, row in enumerate(M)]
    for c in range(n2):
        p = next((r for r in range(c, n2) if A[r][c] != 0), None)
        if p is None:
            return None
        A[c], A[p] = A[p], A[c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        for r in range(n2):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    dens = [A[r][n2 + c].denominator for r in range(n2) for c in range(n2)]
    L = 1
    for d in dens:
        L = L * d // np.gcd(L, d)
    e2 = 0
    while L % (2 ** (e2 + 1)) == 0:
        e2 += 1
    return L, e2

# ---- elements with support <= s and coefficients in {1, i, -1, -i}
UNITS = [1, 1j, -1, -1j]
def elements(maxsupp, coeffs=UNITS):
    out = []
    for s in range(0, maxsupp + 1):
        for supp in itertools.combinations(range(N), s):
            for cs in itertools.product(coeffs, repeat=s):
                a = np.zeros(N, dtype=complex)
                for j, c in zip(supp, cs):
                    a[j] = c
                out.append(a)
    return out

def canon1(a):
    """canonical key of a 1x1 entry up to units (+-1, +-i) x eta^k, Galois eta -> eta^g, complex conjugation"""
    best = None
    for conj in (False, True):
        b0 = np.conj(a) if conj else a
        for g in range(1, N):
            b1 = np.zeros(N, dtype=complex)
            for j in range(N):
                b1[(j * g) % N] = b0[j]
            for k in range(N):
                b2 = np.roll(b1, k)
                for u in UNITS:
                    b3 = b2 * u
                    key = tuple((int(round(x.real)), int(round(x.imag))) for x in b3)
                    if best is None or key < best:
                        best = key
    return best

def fmt(a):
    terms = []
    for j in range(N):
        c = a[j]
        if c == 0:
            continue
        cs = {1: "+", -1: "-", 1j: "+i", -1j: "-i"}.get(complex(c), "+(%s)" % c)
        terms.append("%s%s" % (cs, "1" if j == 0 else ("e" if j == 1 else "e^%d" % j)))
    return "".join(terms) if terms else "0"

def main():
    out = []
    P = lambda *s: out.append(" ".join(str(x) for x in s))
    P("sep_algebra.py: separation entries/matrices over Z[i][Z7], eta of order 7 (Singer cycle of F_2^3 = F_2^4/<x>)")
    P("residue factors of F_2[Z7]: F2 (t=1), F8 (root of t^3+t+1), F8' (root of t^3+t^2+1); a = -1 swaps F8 and F8'")
    # ---------------- 1x1 entries
    P("\n== 1x1 entries: support <= 3, coefficients +-1, +-i ==")
    classes = {}
    stats = Counter()
    for a in elements(3):
        if not a.any():
            continue
        k = canon1(a)
        if k in classes:
            continue
        inv_c = complex_invertible(a)
        pat, f = residue_pattern(a)
        classes[k] = (a, inv_c, pat)
    P("classes up to units, rotation, Galois, conjugation:", len(classes))
    good = []
    for k, (a, inv_c, pat) in sorted(classes.items(), key=lambda t: (np.count_nonzero(t[1][0]), t[0])):
        tag = ("C-invertible" if inv_c else "C-singular") + ", mod-2 kills " + (",".join(pat) if pat else "nothing (unit over Z_(pi))")
        stats[(np.count_nonzero(a), inv_c, pat)] += 1
        if inv_c and pat:
            dd = inverse_denominators([[a]])
            good.append((a, pat, dd))
    P("counts (support, complex-invertible, mod-2 kill pattern) -> number of classes:")
    for key, v in sorted(stats.items(), key=lambda t: (t[0][0], not t[0][1], t[0][2])):
        P("  ", key, v)
    P("\nentries invertible over C AND singular mod 2 (the separation entries), with exact inverse denominators:")
    for a, pat, dd in good:
        P("   %-22s kills %-12s inverse denominator lcm %s (2-adic exponent %s)" % (fmt(a), ",".join(pat), dd[0], dd[1]))
    # Fano line elements: L L^* = 2 + J
    L = np.zeros(N, dtype=complex); L[[0, 1, 3]] = 1
    LL = mul(L, np.conj(np.roll(L[::-1], 1)))
    P("\nFano line indicator L = 1+e+e^3 (lines of PG(2,2) = eta-translates of the difference set {0,1,3}):")
    P("   L * L^dagger =", [int(round(x.real)) for x in LL], "(= 2 + J: |L(chi)|^2 = 2 for chi != 1, L(1) = 3)")
    P("   characters of L:", [round(abs(c), 6) for c in chars(L)], " residue kill:", residue_pattern(L)[0],
      " inverse:", inverse_denominators([[L]]))
    # ---------------- 2x2 matrices
    P("\n== W = 2: matrices [[1, a], [b, c]] (first entry normalised to 1 by row/column units), a, b monomials or 0, c of support <= 2 ==")
    E1 = [e for e in elements(1)]
    E2 = [e for e in elements(2)]
    one = np.zeros(N, dtype=complex); one[0] = 1
    cnt = Counter(); examples = {}
    # det = c - a b ; classify by (support pattern of det, kill pattern)
    tested = 0
    for a in E1:
        for b in E1:
            ab = mul(a, b)
            for c in E2:
                tested += 1
                d = c - ab
                if not d.any():
                    cnt[("C-singular",)] += 1; continue
                inv_c = complex_invertible(d)
                pat = residue_pattern(d)[0]
                key = (inv_c, pat)
                cnt[key] += 1
                sz = sum(np.count_nonzero(x) for x in (a, b, c))
                if inv_c and pat and (key not in examples or sz < examples[key][0]):
                    examples[key] = (sz, a, b, c, d)
    P("matrices tested:", tested)
    for key, v in sorted(cnt.items(), key=lambda t: str(t[0])):
        P("  ", key, v)
    P("smallest example per class (invertible over C, singular mod 2):")
    for key, (sz, a, b, c, d) in sorted(examples.items(), key=lambda t: str(t[0])):
        dd = inverse_denominators([[one, a], [b, c]])
        P("   kills %-14s T = [[1, %s], [%s, %s]]  det = %s  inverse denominators %s" % (",".join(key[1]), fmt(a), fmt(b), fmt(c), fmt(d), dd))
    # ---------------- W = 3, 4: structured families I + e P (P a W-cycle) and Fano circulants
    P("\n== W = 3, 4: structured families ==")
    for W in (2, 3, 4):
        for k in range(1, N):
            for sgn in (1, -1):
                T = [[np.zeros(N, dtype=complex) for _ in range(W)] for _ in range(W)]
                for r in range(W):
                    T[r][r] = one.copy()
                    e = np.zeros(N, dtype=complex); e[k] = sgn
                    T[r][(r + 1) % W] = e
                d = det(T)
                inv_c = complex_invertible(d); pat = residue_pattern(d)[0]
                if k == 1:
                    dd = inverse_denominators(T) if inv_c else None
                    P("  W=%d  T = I %s e^%d P_cycle: det = %-18s C-invertible %s, mod-2 kills %s, inverse dens %s" % (W, "+" if sgn > 0 else "-", k, fmt(d), inv_c, pat, dd))
    # W x W with Fano-line entries: T = I + L * (all-ones off-diagonal)?  record a few
    for W in (2, 3, 4):
        T = [[(L.copy() if r == c else one.copy()) for c in range(W)] for r in range(W)]
        d = det(T); inv_c = complex_invertible(d); pat = residue_pattern(d)[0]
        P("  W=%d  T = L on the diagonal, 1 off it: det = %s, C-invertible %s, kills %s, inverse dens %s" % (W, fmt(d) if np.count_nonzero(d) < 8 else str([complex(x) for x in d]), inv_c, pat, inverse_denominators(T) if inv_c else None))
    # ---------------- RM(1,3) coset matching
    P("\n== matching with the 16 cosets of RM(1,3) = extended Hamming [8,4,4] ==")
    P("syndrome of the 8-bit word (f, parity bit 0) for f = det mod 2 in F_2[Z7]: (f(alpha) in F8, f(1) in F2).")
    P("RM(1,3) cosets <-> (f(alpha), f(1)): the code itself = (0, 0); the 8 weight-1 cosets = (alpha^j, 1) and (0, 1);")
    P("the 7 weight-2 cosets = (alpha^j, 0).  So: kills F8 <-> f(alpha) = 0 <-> coset (0, *); kills F2 <-> (*, 0).")
    tab = Counter()
    for a, pat, dd in good:
        f = mod2(a)
        tab[(ev8(f, 1), sum(f) % 2, ev8(f, 6))] += 1
    P("separation entries (1x1 classes above) by (f(alpha), f(1), f(alpha^-1)):", dict(tab))
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out_sep_algebra.txt")
    open(path, "w").write("\n".join(out) + "\n")
    print("\n".join(out))

if __name__ == "__main__":
    main()
