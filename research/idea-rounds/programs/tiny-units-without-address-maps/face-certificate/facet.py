#!/usr/bin/env python3
"""facet.py m f [mode]
Antisymmetric part A = Lambda^f(C^(2^m)) of f digit columns, basis a(T), T an f-subset of the addresses.
N = span of a(T) over the f-subsets T that lie in a facet {x_j = c} of the m-cube (mode 'facet'),
or, mode 'U:u1,u2,...' : facets {u.x = c} for the listed functionals u (bit masks).
Prints cap_w(N) = rank(dir(w): N -> A/N) for every direction w and dem(N) = rank(kernel: N -> A/N).
Ranks: floating point (singular values) and exactly modulo the prime 998244353 (Gaussian integers, i -> a square root of -1).
"""
import sys, itertools, math
import numpy as np
m = int(sys.argv[1]); f = int(sys.argv[2])
mode = sys.argv[3] if len(sys.argv) > 3 else 'facet'
n = 1 << m
P = 998244353
I_P = pow(3, (P - 1) // 4, P)
assert (I_P * I_P) % P == P - 1

def par(v):
    return bin(v).count('1') & 1

def shift(w):
    X = np.zeros((n, n), dtype=np.complex128)
    for x in range(n):
        X[x ^ w, x] = 1
    return X

def dir2(w):  # 2*dir(w), Gaussian integer entries
    return (1 + 1j) * np.eye(n) + (1 - 1j) * shift(w)

def kernel_scaled():  # 2^m * kernel
    K = np.eye(n, dtype=np.complex128)
    for j in range(m):
        K = dir2(1 << j) @ K
    return K

if mode == 'facet':
    U = [1 << j for j in range(m)]
else:
    U = [int(t) for t in mode.split(':')[1].split(',')]

sets = np.array(list(itertools.combinations(range(n), f)), dtype=np.int64)
def in_N(T):
    for u in U:
        s = sum(par(u & x) for x in T)
        if s == 0 or s == f:
            return True
    return False
flag = np.array([in_N(T) for T in sets])
cols = sets[flag]; rows = sets[~flag]
print(f"m={m} f={f} mode={mode} U={U} dim A={len(sets)} dim N={len(cols)} rows outside N={len(rows)}", flush=True)

def minors(g):
    """matrix of minors det(g[S,T]) for S in rows (outside N), T in cols (N); Gaussian integers as two int64 arrays"""
    R = len(rows); C = len(cols)
    re = np.zeros((R, C), dtype=np.int64); im = np.zeros((R, C), dtype=np.int64)
    worst = 0.0
    for c in range(C):
        sub = g[:, cols[c]]              # n x f
        stack = sub[rows]                # R x f x f
        d = np.linalg.det(stack)
        r_ = np.rint(d.real); i_ = np.rint(d.imag)
        worst = max(worst, float(np.max(np.abs(d.real - r_))), float(np.max(np.abs(d.imag - i_))))
        re[:, c] = r_.astype(np.int64); im[:, c] = i_.astype(np.int64)
    return re, im, worst

def rank_mod_p(re, im):
    M = (re % P + (im % P) * I_P) % P
    M = M.copy()
    R, C = M.shape
    r = 0
    for c in range(C):
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0:
            continue
        p0 = r + piv[0]
        if p0 != r:
            M[[r, p0]] = M[[p0, r]]
        inv = pow(int(M[r, c]), P - 2, P)
        M[r] = (M[r] * inv) % P
        nz = np.nonzero(M[:, c])[0]
        nz = nz[nz != r]
        if len(nz):
            fac = M[nz, c][:, None]
            # products below 2^60
            M[nz] = (M[nz] - (fac * M[r][None, :]) % P) % P
        r += 1
        if r == R:
            break
    return r

def rank_float(re, im):
    M = re.astype(np.float64) + 1j * im.astype(np.float64)
    if min(M.shape) == 0:
        return 0, 0.0, 0.0
    s = np.linalg.svd(M, compute_uv=False)
    tol = s[0] * 1e-9 if s[0] > 0 else 0
    k = int(np.sum(s > tol))
    return k, (float(s[k - 1]) if k else 0.0), (float(s[k]) if k < len(s) else 0.0)

def report(name, g):
    re, im, worst = minors(g)
    rp = rank_mod_p(re, im)
    rf, smin, snext = rank_float(re, im)
    print(f"{name}: rank mod p = {rp}; floating rank = {rf} (least kept singular value {smin:.3g}, next {snext:.3g}); worst rounding of a minor {worst:.2g}", flush=True)
    return rp, rf

done = {}
for w in range(1, n):
    key = bin(w).count('1') if mode == 'facet' else w
    if mode == 'facet' and key in done:
        continue
    done[key] = report(f"cap, direction w={w} (weight {bin(w).count('1')})", dir2(w))
dem = report("dem (kernel)", kernel_scaled())
c = len(cols) // len(U)
print(f"summary m={m} f={f}: dim N={len(cols)} = {len(U)} x {c}; caps={ {k: v[0] for k, v in done.items()} }; dem={dem[0]}; dem/max cap = {dem[0]/max(v[0] for v in done.values()):.4f}")
