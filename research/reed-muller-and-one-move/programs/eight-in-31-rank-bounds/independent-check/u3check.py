# u3check.py (v-r3-classes): my vr3 bounds (DUMP=1 POOL=1, block 3, targets 0,0) against the TRUE E of genuine
# numerical U_3 = 2 in 7 units (ex2-families solver refl2.System2).  Direction z (refl2 coords) <-> my field element
# phi(z) = XOR of f_j over the bits of z, f = (3,5,7) my self-dual basis of GF(8) (mod x^3+x+1).
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, '../../shared')
from refl2 import System2, random_x
f = (3, 5, 7)
def phi(z):
    y = 0
    for j in range(3):
        if z >> j & 1: y ^= f[j]
    return y
inv = {phi(z): z for z in range(1, 8)}
B = set()
for l in open(sys.argv[1]):
    t = l.split()
    if t and t[0] == 'B': B.add((int(t[1]), int(t[2]), int(t[3])))
m, S, W = 3, [1, 2, 3, 4, 5, 6, 7], 2
for seed in (5, 6, 7, 8):
    sy = System2(m, S, W); rng = np.random.default_rng(seed)
    for k in range(8):
        x, fres, it = sy.solve(random_x(sy, rng), iters=600)
        if fres < 1e-20: break
    A_, B_, H = sy.unpack(x); K = B_.T @ A_; E = np.eye(7) + H - K
    bad = tl = tu = 0
    for (mask, lo, hi) in B:
        D = [S.index(inv[y]) for y in range(1, 8) if mask >> y & 1]
        sv = np.linalg.svd(E[np.ix_(D, D)], compute_uv=False)
        rk = int((sv > 1e-7 * max(1, sv[0])).sum())
        if not (lo <= rk <= hi): bad += 1; print("VIOLATION", D, lo, hi, rk)
        tl += rk == lo; tu += rk == hi
    # T10 with my parities
    import struct
    raw = open('tab3.bin', 'rb').read(); nq = struct.unpack('<II', raw[:8])[1]
    t10bad = 0
    for q in range(1, nq):
        mask, par = struct.unpack('<II', raw[8 + 8 * q: 16 + 8 * q])
        r = 2 * (par & 1)
        D = [S.index(inv[y]) for y in range(1, 8) if mask >> y & 1]
        sv = np.linalg.svd(E[np.ix_(D, D)], compute_uv=False)
        rk = int((sv > 1e-7 * max(1, sv[0])).sum())
        t10bad += rk != len(D) - r
    print("seed", seed, "res %.1e" % fres, "distinct bounds", len(B), "violations", bad, "tight l", tl, "tight u", tu, "T10 mismatches", t10bad)
