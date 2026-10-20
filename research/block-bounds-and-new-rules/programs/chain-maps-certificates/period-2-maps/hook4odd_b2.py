# b4-depth4: hook certificate at m = 4 for ODD f (f = k+1, k even), exact mod p, one translation-character block.
# For odd f, Z_c maps the chi-block onto the (chi + c)-block and Z_c K = K (phase X_c Z_c), so
# (N_chi, K N_chi) -> (N_{chi+c}, K N_{chi+c}): all 16 blocks have the same rank; total = 16 x rank(chi = 0).
# The same holds for a single direction D = dir(e1) (Z_c D = D (phase or X-twisted Z_c)), so cap = 16 x cap(chi=0).
# xc-algebra copy of wave16/b4-depth4/hook4odd.py: adds OPS entry 'B2' = dir(1)dir(2) (a rank-2 block); nothing else changed.
# Usage: python3 -I hook4odd.py K  (k even)   ->  cap, dem (mod p rank = lower bound of the true rank)
import sys, os, time, itertools
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np
from math import comb
k = int(sys.argv[1]); m = 4; n = 16
p = 33554393
I = next(pow(g, (p - 1) // 4, p) for g in range(2, 100) if pow(g, (p - 1) // 2, p) == p - 1)
inv2 = (p + 1) // 2; alpha = (1 + I) * inv2 % p; beta = (1 - I) * inv2 % p
def dirmat(z):
    A = np.zeros((n, n), dtype=np.int64)
    for x in range(n):
        A[x, x] = (A[x, x] + alpha) % p; A[x ^ z, x] = (A[x ^ z, x] + beta) % p
    return A
def mm(A, B): return (A @ B) % p
subs = {j: sorted([sum(1 << t for t in c) for c in itertools.combinations(range(n), j)]) for j in range(k + 1)}
idx = {j: {b: i for i, b in enumerate(subs[j])} for j in range(k + 1)}
wtab = {}
for j in range(k):
    Ts = subs[j + 1]
    src = np.zeros((len(Ts), j + 1), dtype=np.int64); sg = np.zeros_like(src); ui = np.zeros_like(src)
    for r, T in enumerate(Ts):
        els = [t for t in range(n) if T >> t & 1]
        for c, i in enumerate(els):
            src[r, c] = idx[j][T & ~(1 << i)]; sg[r, c] = -1 if sum(1 for s in els if s > i) % 2 else 1; ui[r, c] = i
    wtab[j] = (src, sg, ui)
def wedge_cols(A, Smask):
    lam = np.ones(1, dtype=np.int64)
    for j, z in enumerate([t for t in range(n) if Smask >> t & 1]):
        src, sg, ui = wtab[j]; lam = ((lam[src] * A[:, z][ui]) % p * sg).sum(axis=1) % p
    return lam
Sk = subs[k]; NK = len(Sk)
def translate(S, c): return sum(1 << (t ^ c) for t in range(n) if S >> t & 1)
tr = np.array([[idx[k][translate(S, c)] for S in Sk] for c in range(n)], dtype=np.int64)
def sortsign(seq):
    return -1 if sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] > seq[j]) % 2 else 1
eps = np.array([[sortsign([t ^ c for t in range(n) if S >> t & 1]) for S in Sk] for c in range(n)], dtype=np.int64)
inN = np.array([S & 1 == 1 for S in Sk]); Ncols = np.nonzero(inN)[0]; Outrows = np.nonzero(~inN)[0]
print(f"m=4 k={k} f={k+1}: dim W = {n*NK}, dim N = {n*len(Ncols)}, block chi=0: N {len(Ncols)}, out {len(Outrows)}", flush=True)
def block0(A):
    M = np.zeros((len(Outrows), len(Ncols)), dtype=np.int64)
    for j, ci in enumerate(Ncols):
        lam = wedge_cols(A, Sk[ci])
        G = (A[:, 0][:, None] * lam[tr] * eps) % p
        M[:, j] = G.sum(axis=0)[Outrows] % p          # chi = 0: plain sum over c
    return M
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
D = {z: dirmat(z) for z in (1, 2, 4, 8)}
K = mm(mm(mm(D[1], D[2]), D[4]), D[8])
OPS = os.environ.get('OPS', 'd1,K').split(',')
B2 = mm(D[1], D[2])
for name, A in [(nm, {'d1': D[1], 'K': K, 'B2': B2}[nm]) for nm in OPS]:
    t0 = time.time(); M = block0(A); t1 = time.time(); r = rank_modp(M)
    print(f"  {name}: rank on block chi=0 = {r}; x16 = {16*r}  (columns {t1-t0:.0f}s, rank {time.time()-t1:.0f}s)", flush=True)
capf = (n // 2) * comb(n - 2, k - 1); dimN = k * comb(n, k); Hb = n * comb(n, k) - comb(n, k + 1) - dimN
print(f"cap formula {capf}; dim N {dimN}; hook bound (dim H - dim N) {Hb}; ratio bound {min(dimN, Hb)/capf:.4f}")
