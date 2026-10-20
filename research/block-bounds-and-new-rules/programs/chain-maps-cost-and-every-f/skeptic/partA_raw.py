# xc-v-legit (wave 17 skeptic): the conjugated-move identity on RAW ARRAYS, exact (Gaussian integers, int64 with a bound
# check), written from scratch.  Engine layout idx = spec*2^(mf) + sum_b row_b 2^(fb); column j = bit j of every row.
#   C = [[a,b],[b,a]], a=(1+i)/2, b=(1-i)/2 (littleC, BinaryFrames.lean);  D_u = a*1 + b*X_u (X_u: x -> x+u);
#   paid move  M_z = D_z^(x)f = prod_j D_{z e_j}, computed TWO ways: (1) product of f butterfly passes, (2) as ONE
#   "recursive call": gather each coset of span{z e_j} into a length-2^f fibre (coefficient order), apply C^(x)f to the
#   fibre by an independent dense Kronecker product, scatter back  -> checks claim (b) literally.
#   P = row ^= row<<1 on every row as a permutation operator (P g)(Ax) = g(x), A = I + S per row (truncated mod S^f).
# Claim (g-hypotheses 19:16, xc-legit Part A): P M_z P^-1 = X_{z x 1} . D_{z e_0}^-1 . CZ^W_z  for every f, where
#   CZ^W_z = H^(x)k diag((-1)^{sum_j y_j y_{j+1}}) H^(x)k / 2^k, y_j = z . w_j (w = Walsh index, w_j its column j).
# Also: the control (D^-1 on every column), the other orientation (P^-1 M P), and r = 1 spare bits.
# Usage: python3 -I partA_raw.py     (largest: m=2 f=7 r=1 -> 2^15 entries; m=3 f=4 r=0 -> 2^12; < 100 MB)
import numpy as np, itertools, random

BOUND = 1 << 60
class GV:     # Gaussian-integer vector
    def __init__(s, re, im): s.re = np.asarray(re, dtype=np.int64); s.im = np.asarray(im, dtype=np.int64); s.chk()
    def chk(s): assert np.abs(s.re).max(initial=0) < BOUND and np.abs(s.im).max(initial=0) < BOUND
    def perm_from(s, src):            # out[i] = s[src[i]]
        return GV(s.re[src], s.im[src])
    def __add__(s, o): return GV(s.re + o.re, s.im + o.im)
    def __sub__(s, o): return GV(s.re - o.re, s.im - o.im)
    def times(s, c):                  # c = (cr, ci) Gaussian integer
        cr, ci = c; return GV(s.re * cr - s.im * ci, s.re * ci + s.im * cr)
    def scale(s, t): return GV(s.re * t, s.im * t)
    def eq(s, o): return np.array_equal(s.re, o.re) and np.array_equal(s.im, o.im)

def bitpos(m, f):     # address bit (b, j) -> index bit position
    return lambda b, j: f * b + j
def vec_u(m, f, z, j):  # index mask of z (x) e_j
    return sum(1 << (f * b + j) for b in range(m) if (z >> b) & 1)

def D2(g, u, N, inv=False):      # 2*D_u g  (inv: 2*D_u^{-1} g = (1-i) g + (1+i) g(x+u))
    src = np.arange(N) ^ u; gx = g.perm_from(src)
    if not inv: return g.times((1, 1)) + gx.times((1, -1))
    return g.times((1, -1)) + gx.times((1, 1))

def move_butterflies(g, m, f, z, N):     # 2^f * D_z^(x)f g
    for j in range(f): g = D2(g, vec_u(m, f, z, j), N)
    return g

def littleC2():  # 2*C as Gaussian-integer 2x2 (re, im)
    return np.array([[1, 1], [1, 1]]), np.array([[1, -1], [-1, 1]])
def kron_C2(f):  # 2^f * C^(x)f, dense, as (re, im) int64
    Re, Im = np.array([[1]]), np.array([[0]])
    cr, ci = littleC2()
    for _ in range(f):
        Re, Im = np.kron(Re, cr) - np.kron(Im, ci), np.kron(Re, ci) + np.kron(Im, cr)
    return Re, Im

def move_recursive_call(g, m, f, z, N, KC):   # 2^f * D_z^(x)f g via fibres = one size-f call on 2^(k-f) fibres
    us = [vec_u(m, f, z, j) for j in range(f)]
    span = np.zeros(1 << f, dtype=np.int64)
    for c in range(1 << f):
        v = 0
        for j in range(f):
            if (c >> j) & 1: v ^= us[j]
        span[c] = v
    seen = np.zeros(N, dtype=bool); outre = np.zeros(N, np.int64); outim = np.zeros(N, np.int64)
    Kre, Kim = KC
    # C^(x)f index convention: kron builds bit order (c_0 most significant); match by reversing coefficient bits
    rev = np.array([int(format(c, f'0{f}b')[::-1], 2) if f else 0 for c in range(1 << f)])
    fibres = 0
    for base in range(N):
        if seen[base]: continue
        idx = base ^ span; seen[idx] = True; fibres += 1
        xr, xi = g.re[idx][rev], g.im[idx][rev]
        yr = Kre @ xr - Kim @ xi; yi = Kre @ xi + Kim @ xr
        outre[idx[rev]] = yr; outim[idx[rev]] = yi
    return GV(outre, outim), fibres

def chain_perm(m, f, r):     # A: row -> row ^ (row<<1 mod 2^f) on every row; returns array A[x]
    k = m * f + r; N = 1 << k; F = (1 << f) - 1
    x = np.arange(N); out = (x >> (m * f)) << (m * f)
    for b in range(m):
        row = (x >> (f * b)) & F
        out |= (row ^ ((row << 1) & F)) << (f * b)
    return out

def wht(g, N):   # unnormalised Walsh-Hadamard on Gaussian-integer vector
    re, im = g.re.copy(), g.im.copy(); h = 1
    while h < N:
        re = re.reshape(-1, 2, h); im = im.reshape(-1, 2, h)
        re = np.stack([re[:, 0] + re[:, 1], re[:, 0] - re[:, 1]], 1).reshape(-1)
        im = np.stack([im[:, 0] + im[:, 1], im[:, 0] - im[:, 1]], 1).reshape(-1)
        h *= 2
    return GV(re, im)

def czw_signs(m, f, r, z):
    N = 1 << (m * f + r); w = np.arange(N)
    y = np.zeros((f, N), dtype=np.int64)
    for j in range(f):
        acc = np.zeros(N, dtype=np.int64)
        for b in range(m):
            if (z >> b) & 1: acc ^= (w >> (f * b + j)) & 1
        y[j] = acc
    s = np.zeros(N, dtype=np.int64)
    for j in range(f - 1): s ^= y[j] & y[j + 1]
    return 1 - 2 * s

def check(m, f, r, z, rng):
    k = m * f + r; N = 1 << k
    g = GV([rng.randrange(-50, 51) for _ in range(N)], [rng.randrange(-50, 51) for _ in range(N)])
    A = chain_perm(m, f, r)
    assert len(set(A.tolist())) == N            # bijection
    Ainv = np.empty(N, np.int64); Ainv[A] = np.arange(N)
    P = lambda v: v.perm_from(Ainv)              # (P g)(A x) = g(x)  <=> (P g)(y) = g(A^-1 y)
    Pinv = lambda v: v.perm_from(A)
    KC = kron_C2(f)
    # (b): butterflies vs one recursive call per fibre
    mb = move_butterflies(g, m, f, z, N)
    mr, fibres = move_recursive_call(g, m, f, z, N, KC)
    ok_b = mb.eq(mr) and fibres == 1 << (k - f)
    # LHS  2^f * P M P^-1 g
    lhs = P(move_butterflies(Pinv(g), m, f, z, N))
    lhs_other = Pinv(move_butterflies(P(g), m, f, z, N))
    # RHS  X_{z x 1} . 2 D_{z e0}^-1 . 2^k CZ^W   (all three commute; apply right to left)
    s = czw_signs(m, f, r, z)
    t = wht(g, N); t = GV(t.re * s, t.im * s); t = wht(t, N)          # 2^k * CZ^W g
    t = D2(t, vec_u(m, f, z, 0), N, inv=True)                          # * 2
    zall = sum(vec_u(m, f, z, j) for j in range(f))
    rhs = t.perm_from(np.arange(N) ^ zall)                             # X_{z x 1}
    # control: D^-1 on every column instead of column 0
    c = wht(g, N); c = GV(c.re * s, c.im * s); c = wht(c, N)
    for j in range(f): c = D2(c, vec_u(m, f, z, j), N, inv=True)
    ctl = c.perm_from(np.arange(N) ^ zall)
    # compare with matching scales: lhs has 2^f, rhs has 2^(k+1), ctl has 2^(k+f)
    ok_A = lhs.scale(1 << (k + 1)).eq(rhs.scale(1 << f))
    ok_other = lhs_other.scale(1 << (k + 1)).eq(rhs.scale(1 << f))
    ok_ctl = lhs.scale(1 << (k + f)).eq(ctl.scale(1 << f))
    return ok_b, ok_A, ok_other, ok_ctl

rng = random.Random(7)
print("m f r z | move = one size-f call on 2^(k-f) fibres | P M P^-1 == RHS | P^-1 M P == RHS | control(D^-1 every col) == LHS")
allA = True
for (m, fs, rs) in [(2, range(1, 8), (0, 1)), (3, range(1, 5), (0, 1)), (4, range(1, 4), (0,))]:
    for f in fs:
        for r in rs:
            if m * f + r > 15: continue
            for z in range(1, 1 << m):
                ok_b, ok_A, ok_o, ok_c = check(m, f, r, z, rng)
                allA &= ok_b and ok_A
                print(f"{m} {f} {r} {z} | {ok_b} | {ok_A} | {ok_o} | {ok_c}")
print("PART A + (b) HOLD AT ALL TESTED (m,f,r,z)" if allA else "FAILURE")
