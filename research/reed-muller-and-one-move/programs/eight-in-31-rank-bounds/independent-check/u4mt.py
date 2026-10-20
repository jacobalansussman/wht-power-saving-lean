# u4mt.py (v-r3-classes): genuine numerical programs with a TARGET PER ARRAY (T_k = mult by alpha^(2k) in my self-dual
# basis f; frame alpha^k f_j) at block m, all 2^m-1 directions once, found by ex2-families' refl2.System2 with its
# residual changed only in the target term (F - diag(p_l(Q)) instead of F - p(Q) I).  Then every vr3 'B mask l u'
# bound in each given dump is checked against the true E, and T10 with my parities.
# usage: python3 u4mt.py m ks tab.bin restarts iters "seeds" dump1 [dump2 ...]
import sys, struct, itertools, inspect, time, numpy as np, warnings
warnings.filterwarnings("ignore")
import refl2
from refl2 import System2, random_x
from refl import quad_forms
m = int(sys.argv[1]); ks = [int(x) for x in sys.argv[2].split(',')]; tab = sys.argv[3]
restarts, iters = int(sys.argv[4]), int(sys.argv[5]); seeds = [int(s) for s in sys.argv[6].split(',')]
dumps = sys.argv[7:]; W = len(ks)
POLY = {3: 0b1011, 4: 0b10011, 5: 0b100101}[m]; n = 1 << m
def mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & n: a ^= POLY
    return r
def trf(x):
    t, y = 0, x
    for _ in range(m): t ^= y; y = mul(y, y)
    return t
f = next(B for B in itertools.combinations(range(1, n), m)
         if all(trf(mul(B[i], B[j])) == (i == j) for i in range(m) for j in range(m)))
def phi(z):
    y = 0
    for j in range(m):
        if z >> j & 1: y ^= f[j]
    return y
inv = {phi(z): z for z in range(1, n)}
al = [1]
for k in range(n - 2): al.append(mul(al[-1], 2))
S = list(range(1, n))
src = inspect.getsource(System2.residual)
old = "R = F - self.phis[q][:, None, None] * I[None]"
assert src.count(old) == 1
src = src.replace(old, "R = F - self.phisW[q][:, None, :] * I[None]")
src = "\n".join(l[4:] for l in src.splitlines())
ns = {}; exec(src, refl2.__dict__, ns)
class SystemMT(System2):
    residual = ns['residual']
    def __init__(self):
        super().__init__(m, S, W)
        vals, _ = quad_forms(m)
        ph = []
        for q in range(vals.shape[0]):
            if not vals[q, S].any(): continue
            ph.append([sum(vals[q, inv[mul(al[k], fj)]] for fj in f) % 2 for k in ks])
        self.phisW = np.array(ph, dtype=float)
        assert len(ph) == len(self.Ds)
raw = open(tab, 'rb').read(); nq = struct.unpack('<II', raw[:8])[1]
forms = [struct.unpack('<II', raw[8 + 8 * q: 16 + 8 * q]) for q in range(nq)]
def rk(E, mask):
    D = [S.index(inv[y]) for y in range(1, n) if mask >> y & 1]
    if not D: return 0, None
    sv = np.linalg.svd(E[np.ix_(D, D)], compute_uv=False)
    return int((sv > 1e-7 * max(1, sv[0])).sum()), sv
Bsets = {}
for d in dumps:
    B = set()
    for l in open(d):
        t = l.split()
        if t and t[0] == 'B': B.add((int(t[1]), int(t[2]), int(t[3])))
    Bsets[d] = B
print("m", m, "ks", ks, "basis", f, "dumps", {d: len(b) for d, b in Bsets.items()}, flush=True)
for seed in seeds:
    sy = SystemMT(); rng = np.random.default_rng(seed); best = None; t0 = time.time()
    for k in range(restarts):
        x, fres, it = sy.solve(random_x(sy, rng), iters=iters)
        print("  seed", seed, "restart", k, "res %.2e" % fres, "it", it, "t %.0fs" % (time.time() - t0), flush=True)
        if best is None or fres < best[1]: best = (x, fres)
        if fres < 1e-20: break
    x, fres = best
    if fres > 1e-14: print("seed", seed, "NO SOLUTION (best %.2e)" % fres, flush=True); continue
    A_, B_, H = sy.unpack(x); K = B_.T @ A_; E = np.eye(len(S)) + H - K
    t10 = 0
    for q in range(1, nq):
        mask, par = forms[q]; r = sum(par >> k & 1 for k in ks)
        rr, _ = rk(E, mask); t10 += rr != bin(mask).count('1') - r
    for d, Bs in Bsets.items():
        bad = tl = tu = nz = nl = 0
        for (mask, lo, hi) in Bs:
            r, sv = rk(E, mask)
            nz += hi == 0; nl += lo > 0
            if not (lo <= r <= hi): bad += 1; print("VIOLATION", d, "mask", mask, "l", lo, "u", hi, "rank", r, "sv", np.round(sv, 9)[:6])
            tl += r == lo; tu += r == hi
        print("seed", seed, "res %.1e" % fres, d, "bounds", len(Bs), "zero-blocks", nz, "nonzero-l", nl, "violations", bad,
              "tight l", tl, "tight u", tu, "T10 mismatches", t10, flush=True)
