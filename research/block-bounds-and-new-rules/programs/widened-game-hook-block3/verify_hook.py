# d-tof-cert, wave 16.  EXACT check of the hook certificate.
# Module W = V (x) Lambda^(f-1) V, V = C^(2^m) (a GL(V)-submodule of V^(x)f, so a legal Theorem-R module: the engine's
# f-column data restricted to it).  N = span{ e_x (x) e_T : x in T } (a coordinate subspace of W, kept by every monomial
# matrix: all address permutations, all diagonals).  Claims:
#   cap(dir e1) <= c0 := dim N - (n/2) [C(n,f-2) - C(n-2,f-4)]   (PROVED by exhibiting the vectors
#        u = e_a (x) (e_a ^ w) - e_b (x) (e_b ^ w), {a,b} = {a, a+e1}, w in Lambda^(f-2): R u = i * (same form) in N;
#        checked here exactly in Z[i] for every basis w, independence checked by an exact mod-p rank),
#   dem >= rank_p[N | K N] - dim N  (mod-p rank is a lower bound for the rank over Q(i)),
#   all directions are conjugate to dir(e1) by affine address maps; f(rank-2 block) <= 2 cap (triangle inequality).
# Usage: python3 -I verify_hook.py M F
import sys, itertools, time
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np
from math import comb

m, f = int(sys.argv[1]), int(sys.argv[2])
n = 2 ** m; k = f - 1
def is_prime(q):
    if q < 2: return False
    d = 2
    while d * d <= q:
        if q % d == 0: return False
        d += 1
    return True
p = 33554393
while not (is_prime(p) and p % 4 == 1): p -= 2
I = next(pow(g, (p - 1) // 4, p) for g in range(2, 100) if pow(g, (p - 1) // 2, p) == p - 1)
assert I * I % p == p - 1
subsets = list(itertools.combinations(range(n), k)); sidx = {T: i for i, T in enumerate(subsets)}
basis = [(x, T) for x in range(n) for T in subsets]; bidx = {b: i for i, b in enumerate(basis)}
DW = len(basis)
inN = np.array([x in T for (x, T) in basis])
dimN = int(inN.sum())
# exact Gaussian integers as (re, im) python ints; 2*dir(z): e_x -> (1+i) e_x + (1-i) e_{x^z}
A2 = (1, 1); B2 = (1, -1)
def gm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def wedge_apply(z, T):
    """(2 dir z)^{wedge k} e_T as dict sorted-subset -> Gaussian int"""
    out = {}
    for choice in itertools.product((0, 1), repeat=len(T)):
        idxs = [t if c == 0 else t ^ z for t, c in zip(T, choice)]
        if len(set(idxs)) < len(idxs): continue
        coef = (1, 0)
        for c in choice: coef = gm(coef, A2 if c == 0 else B2)
        # sign of sorting permutation
        perm = sorted(range(len(idxs)), key=lambda j: idxs[j])
        sgn = 1
        seen = [False] * len(perm)
        for s in range(len(perm)):
            if seen[s]: continue
            L = 0; j = s
            while not seen[j]: seen[j] = True; j = perm[j]; L += 1
            if L % 2 == 0: sgn = -sgn
        key = tuple(sorted(idxs))
        c0 = out.get(key, (0, 0))
        out[key] = (c0[0] + sgn * coef[0], c0[1] + sgn * coef[1])
    return {kk: v for kk, v in out.items() if v != (0, 0)}
def W_apply(z, vec):
    """(2 dir z) (x) (2 dir z)^{wedge k} applied to a sparse vector {(x,T): gaussian int}"""
    out = {}
    for (x, T), c in vec.items():
        wT = wedge_apply(z, T)
        for (y, cy) in ((x, A2), (x ^ z, B2)):
            for T2, cT in wT.items():
                v = gm(c, gm(cy, cT))
                key = (y, T2); c0 = out.get(key, (0, 0))
                out[key] = (c0[0] + v[0], c0[1] + v[1])
    return {kk: v for kk, v in out.items() if v != (0, 0)}
t0 = time.time()
# ---- exact cap proof: the vectors u(pair, w) ----
pairs = [(a, a ^ 1) for a in range(n) if a < a ^ 1]
U = []
bad = 0
for (a, b) in pairs:
    for w in itertools.combinations(range(n), k - 1):
        u = {}
        for (x, s) in ((a, 1), (b, -1)):
            if x in w: continue
            T = tuple(sorted(w + (x,)))
            # e_x ^ e_w = sign * e_T  with sign of moving x to its sorted place: x placed first then sorted
            idxs = (x,) + w
            inv = sum(1 for j in range(1, len(idxs)) if idxs[j] < x)
            sg = -1 if inv % 2 else 1
            u[(x, T)] = (s * sg, 0)
        if not u: continue
        Ru = W_apply(1, u)
        if any(not (y in T2) for (y, T2) in Ru): bad += 1
        U.append(u)
print(f"m={m} f={f}: dim W = {DW}, dim N = {dimN}; pair vectors u: {len(U)}, with R u outside N: {bad}  ({time.time()-t0:.1f}s)", flush=True)
assert bad == 0
def to_modp(vec):
    v = np.zeros(DW, dtype=np.int64)
    for (key, c) in vec.items():
        v[bidx[key]] = (c[0] + c[1] * I) % p
    return v
def rank_modp(M):
    M = M.copy() % p; r = 0; rows, cols = M.shape
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]
        if piv != r: M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): M[nzr] = (M[nzr] - (col[nzr, None] * M[r][None, :]) % p) % p
        r += 1
    return r
rU = rank_modp(np.array([to_modp(u) for u in U]))
cap_upper = dimN - rU
print(f"independent pair vectors (exact mod-p rank, a lower bound) = {rU}  ->  cap <= {cap_upper}   (formula {dimN - (n//2)*(comb(n,k-1)-comb(n-2,k-3) if k>=3 else dimN - (n//2)*comb(n,k-1))})", flush=True)
# ---- mod-p ranks for cap (lower bound, consistency) and dem (lower bound) ----
Nvecs = [{(x, T): (1, 0)} for (x, T) in basis if x in T]
def apply_seq(zs, vec):
    for z in zs: vec = W_apply(z, vec)
    return vec
Nm = np.array([to_modp(v) for v in Nvecs])
RN = np.array([to_modp(apply_seq([1], v)) for v in Nvecs])
cap_low = rank_modp(np.concatenate([Nm, RN])) - dimN
print(f"cap mod p (lower bound) = {cap_low}  ({time.time()-t0:.1f}s)", flush=True)
KN = np.array([to_modp(apply_seq([1 << j for j in range(m)], v)) for v in Nvecs])
dem_low = rank_modp(np.concatenate([Nm, KN])) - dimN
print(f"dem mod p (lower bound) = {dem_low}  ({time.time()-t0:.1f}s)", flush=True)
from fractions import Fraction
print(f"RESULT m={m} f={f}: cap = {cap_upper} (upper bound proved, lower bound {cap_low}), dem >= {dem_low}: "
      f"every unit needs >= {Fraction(dem_low, cap_upper)} W = {dem_low/cap_upper:.4f} W moves; plain {m} W")
