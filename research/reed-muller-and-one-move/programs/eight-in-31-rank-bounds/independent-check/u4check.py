# u4check.py (v-r3-classes): my vr3 bounds (DUMP=1) against the TRUE E of genuine numerical programs found by
# ex2-families' solver refl2.System2 (targets all I = "0,...,0"), at block m with S = all 2^m - 1 directions.
# usage: python3 u4check.py m W tab.bin dumpfile [seeds...]
# Direction z (refl2 coords) <-> my field element phi(z) = XOR of f_j over the bits of z, f = my self-dual basis
# (mytab.py bases[0]).  Checks: every 'B mask l u' line (l <= rank <= u; u = 0 means a derived zero block), and T10
# (rank E[D_Q] = |D_Q| - r_Q with my parities) on every form.
import sys, struct, itertools, numpy as np, warnings
warnings.filterwarnings("ignore")
from refl2 import System2, random_x
m, W, tab, dump = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
seeds = [int(s) for s in sys.argv[5:]] or [1, 2, 3, 4]
POLY = {3: 0b1011, 4: 0b10011, 5: 0b100101}[m]; n = 1 << m
def mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & n: a ^= POLY
    return r
def tr(x):
    t, y = 0, x
    for _ in range(m): t ^= y; y = mul(y, y)
    return t
f = next(B for B in itertools.combinations(range(1, n), m)
         if all(tr(mul(B[i], B[j])) == (i == j) for i in range(m) for j in range(m)))
print("basis", f)
def phi(z):
    y = 0
    for j in range(m):
        if z >> j & 1: y ^= f[j]
    return y
inv = {phi(z): z for z in range(1, n)}
S = list(range(1, n))
Bs = set()
for l in open(dump):
    t = l.split()
    if t and t[0] == 'B': Bs.add((int(t[1]), int(t[2]), int(t[3])))
raw = open(tab, 'rb').read(); nq = struct.unpack('<II', raw[:8])[1]
forms = [struct.unpack('<II', raw[8 + 8 * q: 16 + 8 * q]) for q in range(nq)]
def rk(E, mask):
    D = [S.index(inv[y]) for y in range(1, n) if mask >> y & 1]
    if not D: return 0
    sv = np.linalg.svd(E[np.ix_(D, D)], compute_uv=False)
    return int((sv > 1e-7 * max(1, sv[0])).sum()), sv
for seed in seeds:
    sy = System2(m, S, W); rng = np.random.default_rng(seed)
    best = None
    for k in range(6):
        x, fres, it = sy.solve(random_x(sy, rng), iters=800)
        if best is None or fres < best[1]: best = (x, fres)
        if fres < 1e-20: break
    x, fres = best
    print("seed", seed, "residual %.2e" % fres, flush=True)
    if fres > 1e-14: continue
    A_, B_, H = sy.unpack(x); K = B_.T @ A_; E = np.eye(len(S)) + H - K
    bad = tl = tu = nz = 0; minsv_bad = []
    for (mask, lo, hi) in Bs:
        r, sv = rk(E, mask)
        if hi == 0: nz += 1
        if not (lo <= r <= hi): bad += 1; print("VIOLATION mask", mask, "l", lo, "u", hi, "rank", r, "sv", np.round(sv, 9)[:6])
        tl += r == lo; tu += r == hi
    t10 = 0
    for q in range(1, nq):
        mask, par = forms[q]; r = W * (par & 1)
        rr, _ = rk(E, mask); t10 += rr != bin(mask).count('1') - r
    print("seed", seed, "bounds", len(Bs), "zero-blocks", nz, "violations", bad, "tight l", tl, "tight u", tu,
          "T10 mismatches", t10, "max|E|", "%.2e" % abs(E).max(), flush=True)
