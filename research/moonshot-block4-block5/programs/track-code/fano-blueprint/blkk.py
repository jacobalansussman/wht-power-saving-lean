# tb4-tight (wave 16): entry tests for the EXACTLY TIGHT block-4 case (2 arrays, blocks of rank 2,2,2,1) at f = 9.
# Module M = V (x) Lambda^8 V, V = F_p^16 (p = 1 mod 4, i = sqrt(-1) mod p), N = span{e_x (x) e_T : x in T}.
# Translation-character block chi = 0 only: for every translation-equivariant A (products of dir z), the block
# matrix M_A[S, T] (S, T = 8-subsets, T containing 0 = N reps, S not containing 0 = 'out' reps) is computed as in
# b4-depth4/hook4odd.py (copied logic, own file).  Z_c (diag (-1)^{c.x}) is free, maps the chi-block onto the
# (chi+c)-block, and keeps E = N cap K^-1 N and N_D = N cap D^-1 N for every D = product of dirs (K Z_c K^-1 and
# D Z_c D^-1 are Paulis times scalars), so "E inside N_D" holds iff it holds on the chi = 0 block [proved in NOTE].
# Ranks mod p are lower bounds of the true ranks; a containment found mod p is checked at a second prime.
# Usage: python3 -I blk9.py E  r            -> E_r = null(M_{K_r}) on chi=0, saved as E_r_p<p>.npy
#        python3 -I blk9.py test r D1,D2,..  -> dim(E_r cap N_D) for each named D (see DNAMES)
#        python3 -I blk9.py ptest r D1,..    -> P-side: rank [M_{K_r^-1} | M_{D^-1}] vs rank M_{K_r^-1}
import sys, os, time, itertools
import numpy as np
k = int(os.environ.get('HOOKK', '8')); n = 16
P = int(os.environ.get('PRIME', '33554393'))
p = P
I = next(pow(g, (p - 1) // 4, p) for g in range(2, 100) if pow(g, (p - 1) // 2, p) == p - 1)
inv2 = (p + 1) // 2; alpha = (1 + I) * inv2 % p; beta = (1 - I) * inv2 % p
def dirmat(z, inv=False):
    a, b = (beta, alpha) if inv else (alpha, beta)      # dir^-1 = conj(dir) = beta I + alpha S_z
    A = np.zeros((n, n), dtype=np.int64)
    for x in range(n):
        A[x, x] = (A[x, x] + a) % p; A[x ^ z, x] = (A[x ^ z, x] + b) % p
    return A
def mm(A, B): return (A @ B) % p
def prod(zs, inv=False):
    A = np.eye(n, dtype=np.int64)
    for z in zs: A = mm(A, dirmat(z, inv))
    return A
# K_r: kernel on a basis a1..a4 with e1 = a1 + ... + a_r  (r = 1: the standard kernel)
KBASES = {1: (1, 2, 4, 8), 2: (3, 2, 4, 8), 3: (7, 2, 4, 8), 4: (15, 2, 4, 8)}
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
def block0(A, rows=None, cols=None):
    """chi = 0 block of rho(A): columns = N reps (default) or out reps, rows = out reps (default) or N reps,
    for translation-equivariant A"""
    rows = Outrows if rows is None else rows; cols = Ncols if cols is None else cols
    M = np.zeros((len(rows), len(cols)), dtype=np.int64)
    a0 = A[:, 0]
    for j, ci in enumerate(cols):
        lam = wedge_cols(A, Sk[ci])
        G = (a0[:, None] * lam[tr] * eps) % p
        M[:, j] = G.sum(axis=0)[rows] % p
    return M
def rref(M):
    """reduced row echelon form mod p (in place copy); returns (R, pivots)"""
    M = M.copy() % p; r = 0; rows, cols = M.shape; piv = []
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        pr = r + nz[0]
        if pr != r: M[[r, pr]] = M[[pr, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): M[nzr] = (M[nzr] - (col[nzr, None] * M[r][None, :]) % p) % p
        piv.append(c); r += 1
    return M[:r], piv
def rank(M): return len(rref(M)[1])
def nullspace(M):
    R, piv = rref(M); cols = M.shape[1]
    free = [c for c in range(cols) if c not in set(piv)]
    B = np.zeros((len(free), cols), dtype=np.int64)
    for i, fc in enumerate(free):
        B[i, fc] = 1
        for r, pc in enumerate(piv): B[i, pc] = (-R[r, fc]) % p
    return B
def matmul_modp(A, B):
    # A (r x s) int64 < p, B (s x t); p < 2^25 and s <= 6435 keeps every partial sum below 2^63
    assert p < (1 << 25) and A.shape[1] <= 7000
    return (A @ B) % p
DNAMES = {   # name -> list of directions (product, applied right to left does not matter: they commute)
    'd1': [1], 'd3': [3], 'd7': [7], 'd15': [15],
    'B12': [1, 2], 'B13': [1, 3], 'B16': [1, 6], 'B1_14': [1, 14], 'B3_12': [3, 12], 'B3_5': [3, 5], 'B7_8': [7, 8],
    'B7_9': [7, 9], 'B15_1': [15, 1], 'B3_6': [3, 6],
}
if __name__ == '__main__':
    mode = sys.argv[1]; r = int(sys.argv[2]); t0 = time.time()
    tag = f"r{r}_p{p}"
    if mode == 'E':
        MK = block0(prod(KBASES[r])); t1 = time.time()
        E = nullspace(MK)
        print(f"K_{r} basis {KBASES[r]}: rank M_K = {MK.shape[1]-len(E)}, dim E_0 = {len(E)}  (cols {t1-t0:.0f}s, null {time.time()-t1:.0f}s)", flush=True)
        chk = matmul_modp(MK, E.T); assert not chk.any()
        np.save(f"E_{tag}.npy", E)
    elif mode == 'EP':
        # E_P (+) L_0 = null of the N-rows x out-cols block of K (L is G-invariant with zero N-part)
        MK = block0(prod(KBASES[r]), rows=Ncols, cols=Outrows); t1 = time.time()
        EP = nullspace(MK)
        print(f"K_{r}: rank of N x out block = {MK.shape[1]-len(EP)}, dim (E_P + L)_0 = {len(EP)} (L_0 = 715)  (cols {t1-t0:.0f}s, null {time.time()-t1:.0f}s)", flush=True)
        np.save(f"EP_{tag}.npy", EP)
    elif mode == 'test':
        E = np.load(f"E_{tag}.npy")
        for name in sys.argv[3].split(','):
            t1 = time.time(); MD = block0(prod(DNAMES[name])); rD = rank(MD)
            X = matmul_modp(MD, E.T); rx = rank(X.T)
            print(f"K_{r}, D={name} {DNAMES[name]}: rank M_D = {rD} (x16 = {16*rD}); rank(M_D E_0) = {rx} -> dim(E_0 cap N_D) = {len(E)-rx} of {len(E)}; E inside N_D: {rx == 0}  ({time.time()-t1:.0f}s)", flush=True)
    elif mode == 'ptest':
        MKi = block0(prod(KBASES[r], inv=True)); rK = rank(MKi)
        print(f"K_{r}^-1: rank {rK}", flush=True)
        for name in sys.argv[3].split(','):
            t1 = time.time(); MDi = block0(prod(DNAMES[name], inv=True)); rD = rank(MDi)
            rj = rank(np.concatenate([MKi, MDi], axis=1).T)   # rank of the column span = rank of the transpose
            print(f"K_{r}, D={name}: rank M_D^-1 = {rD}; rank [M_K^-1 | M_D^-1] = {rj}; D^-1 N inside N + K^-1 N: {rj == rK}; excess {rj-rK}  ({time.time()-t1:.0f}s)", flush=True)
