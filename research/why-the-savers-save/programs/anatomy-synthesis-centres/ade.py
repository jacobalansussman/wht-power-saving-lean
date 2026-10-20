"""an-synth calc/ade.py -- the centre of every known saver as an index code on a 'two-field design'.

Our own code, stdlib only, exact integer arithmetic. Run: python3 -I calc/ade.py

A two-field design (OpenAI's Alon triples, our E8, the hub's P112) is a set of n ports with
  * real (or complex) unit vectors f_S: the scatter-gather product A = G/|f|^2 (A_SS = 1);
  * binary labels t_S in F_2^k of odd weight, with t_S . t_T = 1 (mod 2)  ==>  A_ST = 0.
Pairs with A_ST != 0 are 'junk' pairs; they must be F_2-orthogonal so that their carriers only climb.
The number of totals (centres) is t = rank A; the copy charge is about t * k (totals at the full window).

For a simply-laced root system the rule 'B_ST = 1 iff <r_S, r_T> is even' (B = J + G mod 2) gives such
labels automatically: B is symmetric with unit diagonal, so B = L L^T over F_2 with k = rank_F2(B)
columns, and B_ST = 0 exactly on the pairs with <r_S, r_T> = +-1, i.e. exactly where A_ST = +-1/2.

This script, for A_n, D_n, E6, E7, E8 and for the 60 complex lines of E8 viewed as a Z[i]-lattice
(the reflection lines of the complex group G31 in C^4):
  n lines, t = rank, junk per line s, tightness (A^2 = (n/t) A), window k = rank_F2(B0), copy ratio t*k/n.
"""
from itertools import combinations, product
from fractions import Fraction
import json, sys


def rank_q(rows):
    """rank over Q of an integer matrix (list of lists)."""
    M = [[Fraction(x) for x in r] for r in rows]
    rk, ncol = 0, len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def rank_f2(bitrows):
    """rank over F_2 of rows given as Python ints."""
    basis = {}
    for r in bitrows:
        while r:
            h = r.bit_length() - 1
            if h in basis:
                r ^= basis[h]
            else:
                basis[h] = r
                break
    return len(basis)


def canon(v):
    for x in v:
        if x != 0:
            return tuple(v) if x > 0 else tuple(-y for y in v)
    raise ValueError


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def lines_A(n):
    d = n + 1
    out = []
    for i, j in combinations(range(d), 2):
        v = [0] * d; v[i] = 1; v[j] = -1
        out.append(tuple(v))
    return out, 2  # vectors, norm^2


def lines_D(n):
    out = []
    for i, j in combinations(range(n), 2):
        for s in (1, -1):
            v = [0] * n; v[i] = 1; v[j] = s
            out.append(tuple(v))
    return out, 2


def roots_E8_scaled():
    """E8 roots scaled by 2 (norm^2 8): D8 part (+-2,+-2,0..) and (+-1)^8 with an even number of minus signs."""
    R = set()
    for i, j in combinations(range(8), 2):
        for si, sj in product((2, -2), repeat=2):
            v = [0] * 8; v[i] = si; v[j] = sj
            R.add(tuple(v))
    for signs in product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            R.add(tuple(signs))
    assert len(R) == 240
    return sorted(R)


def lines_from_roots(roots):
    return sorted({canon(r) for r in roots})


def analyse(name, vecs, norm2):
    n = len(vecs)
    G = [[dot(a, b) for b in vecs] for a in vecs]
    t = rank_q([list(v) for v in vecs])
    # junk: |<r,s>| = norm2/2 (normalised +-1/2); orthogonal: 0
    half = norm2 // 2
    s_rows = [sum(1 for j in range(n) if j != i and abs(G[i][j]) == half) for i in range(n)]
    other = [sum(1 for j in range(n) if j != i and G[i][j] not in (0, half, -half)) for i in range(n)]
    assert max(other) == 0, name
    # tightness: A = G/norm2; A^2 = c A  <=>  G^2 = c*norm2*G
    G2 = [[sum(G[i][k] * G[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    c = Fraction(G2[0][0], G[0][0] * norm2)
    tight = all(G2[i][j] == c * norm2 * G[i][j] for i in range(n) for j in range(n))
    # binary window: B0_ST = 1 iff <r_S,r_T> even in units of norm2/2 (i.e. A_ST = 0 or diagonal)
    rows = []
    for i in range(n):
        r = 0
        for j in range(n):
            if (G[i][j] // half) % 2 == 0:
                r |= 1 << j
        rows.append(r)
    k = rank_f2(rows)
    s = sum(s_rows) / n
    return dict(design=name, n=n, t=t, junk_per_line=s, junk_min=min(s_rows), junk_max=max(s_rows),
                n_over_t=float(Fraction(n, t)), one_plus_s_over_4=1 + s / 4, tight=tight,
                tight_const=float(c), window_k=k, copy_ratio=t * k / n)


def complex_E8():
    """E8 as a Z[i]-lattice: z_k = x_{2k-1} + i x_{2k}; multiplication by i is J(x1,x2,..) = (-x2,x1,..).
    J preserves E8 (each coordinate pair flips the parity of minus signs; four pairs keep it).
    Returns the 60 complex lines as Gaussian-integer 4-vectors (scaled coordinates)."""
    roots = roots_E8_scaled()
    rs = set(roots)
    J = lambda v: tuple(c for k in range(4) for c in (-v[2 * k + 1], v[2 * k]))
    assert all(J(r) in rs for r in roots)
    seen, lines = set(), []
    for r in roots:
        if r in seen:
            continue
        orb = {r, J(r), J(J(r)), J(J(J(r)))}
        seen |= orb
        lines.append(min(orb))
    assert len(lines) == 60
    # Gaussian vectors: z_k = (re, im)
    return [tuple((v[2 * k], v[2 * k + 1]) for k in range(4)) for v in lines]


def herm(z, w):
    """<z,w> = sum conj(z_k) w_k for Gaussian-integer vectors; returns (re, im)."""
    re = sum(a * c + b * d for (a, b), (c, d) in zip(z, w))
    im = sum(a * d - b * c for (a, b), (c, d) in zip(z, w))
    return re, im


def analyse_complex():
    Z = complex_E8()
    n = len(Z)
    H = [[herm(a, b) for b in Z] for a in Z]
    norm2 = H[0][0][0]  # 8 in scaled coordinates
    assert all(H[i][i] == (norm2, 0) for i in range(n))
    mag = [[H[i][j][0] ** 2 + H[i][j][1] ** 2 for j in range(n)] for i in range(n)]
    # normalised |A_ST|^2 = mag / norm2^2  in {0, 1/4, 1/2, 1}
    hist = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                q = Fraction(mag[i][j], norm2 ** 2)
                hist[q] = hist.get(q, 0) + 1
    hist = {str(k): v // n for k, v in sorted(hist.items())}
    # complex rank = (real rank of realification)/2
    real_rows = []
    for z in Z:
        real_rows.append([c for (a, b) in z for c in (a, b)])
        real_rows.append([c for (a, b) in z for c in (-b, a)])  # i*z
    t = rank_q(real_rows) // 2
    # tightness: A^2 = c A with A = H / norm2, Gaussian matrix product
    def mul(x, y):
        return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
    tight = True
    cconst = None
    for i in range(n):
        for j in range(n):
            acc = [0, 0]
            for k in range(n):
                p = mul(H[i][k], H[k][j]); acc[0] += p[0]; acc[1] += p[1]
            # compare acc with c*norm2*H[i][j]
            if H[i][j] == (0, 0):
                if acc != [0, 0]:
                    tight = False
            else:
                hr, hi = H[i][j]
                # c = acc / (norm2 * H) must be a real constant
                # acc * conj(H) / (norm2 * |H|^2)
                num_re = acc[0] * hr + acc[1] * hi
                num_im = acc[1] * hr - acc[0] * hi
                if num_im != 0:
                    tight = False
                cc = Fraction(num_re, norm2 * (hr * hr + hi * hi))
                if cconst is None:
                    cconst = cc
                elif cc != cconst:
                    tight = False
    rows = []
    for i in range(n):
        r = 0
        for j in range(n):
            if H[i][j] == (0, 0) or i == j:
                r |= 1 << j
        rows.append(r)
    k = rank_f2(rows)
    junk = [sum(1 for j in range(n) if j != i and H[i][j] != (0, 0)) for i in range(n)]
    s = sum(junk) / n
    sumsq = Fraction(sum(mag[0][j] for j in range(n) if j != 0), norm2 ** 2)
    return dict(design="E8 as 60 complex lines in C^4 (G31 lines, Gaussian)", n=n, t=t,
                junk_per_line=s, abs2_histogram_per_line=hist, one_plus_sum_abs2=float(1 + sumsq),
                n_over_t=n / t, tight=tight, tight_const=float(cconst), window_k=k, copy_ratio=t * k / n,
                labels_needed_min_k=(n - 1).bit_length() + 1)


def main():
    res = []
    for nn in (3, 4, 5, 6, 7, 8, 9, 12, 16):
        v, nm = lines_A(nn); res.append(analyse(f"A{nn}", v, nm))
    for nn in (4, 5, 6, 8, 9, 12, 16, 17):
        v, nm = lines_D(nn); res.append(analyse(f"D{nn}", v, nm))
    E8 = roots_E8_scaled()
    r0 = E8[0]
    E7 = [r for r in E8 if dot(r, r0) == 0]
    r1 = next(r for r in E8 if dot(r, r0) == -4)
    E6 = [r for r in E7 if dot(r, r1) == 0]
    assert len(E7) == 126 and len(E6) == 72
    for name, R in (("E6", E6), ("E7", E7), ("E8", E8)):
        res.append(analyse(name, lines_from_roots(R), 8))
    res.append(analyse_complex())
    hdr = f"{'design':52s} {'n':>4s} {'t':>3s} {'junk/line':>9s} {'n/t':>7s} {'1+s/4':>7s} {'tight':>5s} {'k':>3s} {'t*k/n':>7s}"
    print(hdr)
    for r in res:
        print(f"{r['design']:52s} {r['n']:4d} {r['t']:3d} {r['junk_per_line']:9.2f} {r['n_over_t']:7.2f} "
              f"{r.get('one_plus_s_over_4', r.get('one_plus_sum_abs2')):7.2f} {str(r['tight']):>5s} "
              f"{r['window_k']:3d} {r['copy_ratio']:7.3f}")
    cx = res[-1]
    print("complex E8: |A_ST|^2 histogram per line:", cx['abs2_histogram_per_line'],
          " minimum window for 60 distinct odd-weight labels:", cx['labels_needed_min_k'])
    json.dump(res, open(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", "w"), indent=1)


if __name__ == "__main__":
    main()
