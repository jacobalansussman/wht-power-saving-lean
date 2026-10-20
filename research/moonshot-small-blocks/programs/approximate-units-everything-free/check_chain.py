"""v-w3-border: independent checks (own code, no xw2lib).
Address at block 2: column j has bits (b0, b1) = bits 2j, 2j+1 of x; f columns; N = 4^f.
Paid D_z = dir(z)^{(x)f}, dir(z) = a I + b X^z, a = (1+i)/2, b = (1-i)/2.  We use the scaled 2^f D_z (Gaussian integers),
vectors as exact Python-int-valued numpy object arrays? -> use int64 real/imag pairs (values stay < 2^40 here).
C1 chain map a1 (bit 1 of column j ^= bit 1 of column j-1) commutes with D1 (direction e1 = bit 0), at f = 1..8.
   a0 commutes with D2.  Negative controls: a0 vs D1, a1 vs D2, a1 vs D3.
C2 order of a1 grows with f (so a1 generates an infinite group of families; not in any finite H).
C3 a1 is not a product of per-column (or first/last-column) letters.
C4 R7 (w3-border S4) rebuilt: t = pi+ x, u_a = D2 t, u = D1(t + lZ1 u_a), u_b = D2(pi+ u), Psi = pi- u - i lZ1 lX1 u_b = 0.
"""
import numpy as np, sys
rng = np.random.default_rng(7)

def cols(f):
    x = np.arange(4 ** f, dtype=np.int64)
    return x

def perm_a(f, r):
    """running XOR on row r: bit r of column j ^= bit r of column j-1 (j = 1..f-1, original values)"""
    x = cols(f); y = x.copy()
    for j in range(1, f):
        src = (x >> (2 * (j - 1) + r)) & 1
        y ^= src << (2 * j + r)
    return y          # (P v)[y[x]] = v[x]

def apply_perm(sig, v):
    out = np.empty_like(v); out[sig] = v; return out

class G:  # gaussian-int vector: (re, im) int64 arrays (python ints via object for safety)
    pass

def gmul_scalar(re, im, cr, ci):
    return re * cr - im * ci, re * ci + im * cr

def apply_D2x(z, re, im, f):
    """scaled: per column (1+i) v + (1-i) v[x ^ z_col]"""
    x = cols(f)
    for j in range(f):
        sh = x ^ (z << (2 * j))
        r1, i1 = gmul_scalar(re, im, 1, 1)
        r2, i2 = gmul_scalar(re[sh], im[sh], 1, -1)
        re, im = r1 + r2, i1 + i2
    return re, im

def randvec(f):
    n = 4 ** f
    return rng.integers(-3, 4, n).astype(object), rng.integers(-3, 4, n).astype(object)

def eq(a, b):
    return all(np.array_equal(p, q) for p, q in zip(a, b))

def commute_test(f, sig, z, nvec=3):
    ok = True
    for _ in range(nvec):
        re, im = randvec(f)
        A = apply_D2x(z, apply_perm(sig, re), apply_perm(sig, im), f)
        r2, i2 = apply_D2x(z, re, im, f)
        B = (apply_perm(sig, r2), apply_perm(sig, i2))
        ok &= eq(A, B)
    return ok

print("C1: commutation of chain maps with paid letters (exact Gaussian-integer arithmetic, 3 random vectors)")
for f in range(1, 9):
    a0, a1 = perm_a(f, 0), perm_a(f, 1)
    print(f" f={f}: a1D1=D1a1 {commute_test(f,a1,1)}  a0D2=D2a0 {commute_test(f,a0,2)}"
          f"  | controls (expect False for f>=2): a0D1 {commute_test(f,a0,1,1)} a1D2 {commute_test(f,a1,2,1)} a1D3 {commute_test(f,a1,3,1)}")
    sys.stdout.flush()

print("C1b: a1 D1 = D1 a1 on ALL basis vectors (dense check), f = 1..4")
for f in range(1, 5):
    n = 4 ** f; a1 = perm_a(f, 1); good = True
    for k in range(n):
        re = np.zeros(n, dtype=object); im = np.zeros(n, dtype=object); re[k] = 1
        A = apply_D2x(1, apply_perm(a1, re), apply_perm(a1, im), f)
        r2, i2 = apply_D2x(1, re, im, f)
        good &= eq(A, (apply_perm(a1, r2), apply_perm(a1, i2)))
    print(f" f={f}: {good}")

print("C2: order of a1 as a permutation, f = 1..8")
for f in range(1, 9):
    s = perm_a(f, 1); p = s.copy(); k = 1
    ident = np.arange(4 ** f)
    while not np.array_equal(p, ident):
        p = s[p]; k += 1
    print(f" f={f}: order {k}")

print("C3: a1 is not a product of single-column letters: bit 1 of column 1 depends on column 0 (f = 3):",
      any(((perm_a(3, 1)[x] >> 3) & 1) != ((x >> 3) & 1) for x in range(64)))

# ---- C4: R7 rebuilt with own letters (scaled exactly) ----
def X1(v, f):     # X on bit 0 of every column
    x = cols(f); m = 0
    for j in range(f): m |= 1 << (2 * j)
    return v[x ^ m]

def apply_lZ1(re, im, f):  # Z on bit 0 of last column
    x = cols(f); s = 1 - 2 * ((x >> (2 * (f - 1))) & 1)
    return re * s, im * s

def apply_lX1(re, im, f):
    x = cols(f); sh = x ^ (1 << (2 * (f - 1)))
    return re[sh], im[sh]

def pi(sign, re, im, f):   # 2*pi^{+-} = 1 +- X1 (scaled by 2)
    return re + sign * X1(re, f), im + sign * X1(im, f)

print("C4: R7 rebuilt (scaled by 2^(3f) * 8 to stay integral), f = 1..6")
for f in range(1, 7):
    worst = 0; nz = 0
    for _ in range(3):
        re, im = randvec(f)
        tr, ti = pi(+1, re, im, f)                         # 2 t
        uar, uai = apply_D2x(2, tr, ti, f)                  # 2^f * 2 * u_a
        lr, li = apply_lZ1(uar, uai, f)                     # lZ1 u_a (scaled 2^(f+1))
        cr, ci = tr * 2 ** f + lr, ti * 2 ** f + li          # 2^(f+1)(t + lZ1 u_a)
        ur, ui = apply_D2x(1, cr, ci, f)                    # 2^(2f+1) u
        pr, pim = pi(+1, ur, ui, f)                         # 2^(2f+2) pi+ u
        ubr, ubi = apply_D2x(2, pr, pim, f)                 # 2^(3f+2) u_b
        mr, mi = pi(-1, ur, ui, f)                          # 2^(2f+2) pi- u
        # -i lZ1 lX1 u_b: lX1 first then lZ1 (operator product lZ1 lX1 acts lX1 first)
        xr, xi = apply_lX1(ubr, ubi, f); zr, zi = apply_lZ1(xr, xi, f)
        # multiply by -i: (zr + i zi)(-i) = zi - i zr
        tr2, ti2 = zi, -zr
        Pr = mr * 2 ** f + tr2; Pi = mi * 2 ** f + ti2       # 2^(3f+2) Psi
        worst = max(worst, max(abs(int(v)) for v in np.concatenate([Pr, Pi])))
        nz = max(nz, max(abs(int(v)) for v in np.concatenate([mr, mi])))
    print(f" f={f}: max|Psi| (scaled) = {worst}; max|pi- u| (scaled) = {nz} (nonzero half)")
