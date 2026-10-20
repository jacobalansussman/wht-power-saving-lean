# v28-theorem (skeptic), wave 16, 13 Oct 2026.  INDEPENDENT re-computation of the hook certificate of d-tof-cert.
# Written from scratch (no code of d-tof-cert imported).  Module Wm = V (x) Lambda^k V, V = C^(2^m), k = f-1,
# g acts as g (x) Lambda^k g, with Lambda^k g [S,T] = det g[S,T] (rows S, columns T).
# N = span{e_x (x) e_T : x in T}, a coordinate subspace.  For any g:  f(g) = dim(N + gN) - dim N = rank g_W[notN, N].
# Ranks: (a) modulo two primes p = 1 mod 4 (i -> sqrt(-1) mod p): LOWER bounds for the rank over Q(i);
#        (b) EXACT over Z[i] by fraction-free (Bareiss) elimination with Python integers: the true rank.
# Matrices used (scaled to Z[i]; scaling does not change spans):
#   2 dir z = (1+i) I + (1-i) shift z ;  block of rank 2 = (2 dir z1)(2 dir z2) ;  kernel = prod_j (2 dir e_j).
# Usage: python3 -I hookcheck.py M F [exact:0/1] [allblocks:0/1]
import sys, itertools, time

m, f = int(sys.argv[1]), int(sys.argv[2])
EXACT = int(sys.argv[3]) if len(sys.argv) > 3 else 1
ALLB = int(sys.argv[4]) if len(sys.argv) > 4 else 0
n = 1 << m; k = f - 1
t0 = time.time()

# ---------- Gaussian integers as (re, im) ----------
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gadd(a, b): return (a[0]+b[0], a[1]+b[1])
def gsub(a, b): return (a[0]-b[0], a[1]-b[1])
Z0 = (0, 0); ONE = (1, 0)

def matmul(A, B):
    nn = len(A)
    return [[ (sum(gmul(A[i][l], B[l][j])[0] for l in range(nn)), sum(gmul(A[i][l], B[l][j])[1] for l in range(nn)))
              for j in range(nn)] for i in range(nn)]
def dir2(z):
    M = [[Z0]*n for _ in range(n)]
    for x in range(n):
        M[x][x] = gadd(M[x][x], (1, 1))
        M[x][x ^ z] = gadd(M[x][x ^ z], (1, -1))   # (shift z v)(x) = v(x+z): row x, column x^z
    return M

def det(M):
    # Gaussian-integer determinant by Laplace expansion (size <= 5 here)
    s = len(M)
    if s == 0: return ONE
    if s == 1: return M[0][0]
    tot = Z0
    for j in range(s):
        if M[0][j] == Z0: continue
        minor = [row[:j] + row[j+1:] for row in M[1:]]
        d = gmul(M[0][j], det(minor))
        tot = gadd(tot, d) if j % 2 == 0 else gsub(tot, d)
    return tot

subsets = list(itertools.combinations(range(n), k))
basis = [(x, T) for x in range(n) for T in subsets]
inN = [x in T for (x, T) in basis]
Nidx = [i for i, b in enumerate(basis) if inN[i]]
Cidx = [i for i, b in enumerate(basis) if not inN[i]]
dimW, dimN = len(basis), len(Nidx)
sidx = {T: i for i, T in enumerate(subsets)}

def wedge(g):
    """Lambda^k g as dict (S index, T index) -> det g[S,T]  (only nonzero)"""
    out = {}
    for T in subsets:
        # rows that can be nonzero: union of supports of columns in T
        rowsupp = sorted(set(r for t in T for r in range(n) if g[r][t] != Z0))
        for S in itertools.combinations(rowsupp, k):
            d = det([[g[s][t] for t in T] for s in S])
            if d != Z0: out[(sidx[S], sidx[T])] = d
    return out

def block_matrix_cols(g):
    """columns of g_W restricted to N columns, rows restricted to not-N: returns list (over N columns) of dict row->val"""
    lw = wedge(g)
    byT = {}
    for (si, ti), d in lw.items(): byT.setdefault(ti, []).append((si, d))
    rowpos = {ci: r for r, ci in enumerate(Cidx)}
    cols = []
    for j in Nidx:
        x, T = basis[j]
        col = {}
        for y in range(n):
            gy = g[y][x]
            if gy == Z0: continue
            for si, d in byT.get(sidx[T], []):
                S = subsets[si]
                if y in S: continue          # lands in N: projected away
                bi = y * len(subsets) + si
                assert basis[bi] == (y, S)
                r = rowpos[bi]
                v = gmul(gy, d)
                col[r] = gadd(col.get(r, Z0), v)
        cols.append({r: v for r, v in col.items() if v != Z0})
    return cols

def is_prime(q):
    if q < 2: return False
    d = 2
    while d * d <= q:
        if q % d == 0: return False
        d += 1
    return True
def primes_1mod4(start, cnt):
    out = []; p = start
    while len(out) < cnt:
        if is_prime(p) and p % 4 == 1: out.append(p)
        p -= 1
    return out

def rank_modp(cols, nrows, p):
    I = next(pow(g, (p - 1) // 4, p) for g in range(2, 200) if pow(g, (p - 1) // 2, p) == p - 1)
    assert I * I % p == p - 1
    # rows = columns of the matrix (rank is the same); dense list-of-lists mod p
    rows = []
    for col in cols:
        v = [0] * nrows
        for r, (a, b) in col.items(): v[r] = (a + b * I) % p
        rows.append(v)
    rank = 0; ncols = nrows
    piv_rows = []
    for c in range(ncols):
        pr = None
        for i in range(rank, len(rows)):
            if rows[i][c] % p: pr = i; break
        if pr is None: continue
        rows[rank], rows[pr] = rows[pr], rows[rank]
        inv = pow(rows[rank][c], p - 2, p)
        rr = [(v * inv) % p for v in rows[rank]]; rows[rank] = rr
        for i in range(rank + 1, len(rows)):
            a = rows[i][c]
            if a:
                ri = rows[i]
                rows[i] = [(ri[t] - a * rr[t]) % p for t in range(ncols)]
        rank += 1
    return rank

def gdiv_exact(a, b):
    nb = b[0]*b[0] + b[1]*b[1]
    num = gmul(a, (b[0], -b[1]))
    assert num[0] % nb == 0 and num[1] % nb == 0, "Bareiss division not exact"
    return (num[0] // nb, num[1] // nb)

def rank_exact(cols, nrows):
    """fraction-free Gaussian elimination over Z[i] (Bareiss).  Rows of the work matrix = the given columns."""
    M = []
    for col in cols:
        v = [Z0] * nrows
        for r, val in col.items(): v[r] = val
        M.append(v)
    R = len(M); rank = 0; prev = ONE
    for c in range(nrows):
        pr = None
        for i in range(rank, R):
            if M[i][c] != Z0: pr = i; break
        if pr is None: continue
        M[rank], M[pr] = M[pr], M[rank]
        piv = M[rank][c]; prow = M[rank]
        for i in range(rank + 1, R):
            row = M[i]; a = row[c]
            new = [Z0] * nrows
            for t in range(c, nrows):
                v = gsub(gmul(piv, row[t]), gmul(a, prow[t]))
                new[t] = gdiv_exact(v, prev) if v != Z0 else Z0
            M[i] = new
        prev = piv; rank += 1
    return rank

P = primes_1mod4(2**31 - 1, 2)
def fval(name, g):
    cols = block_matrix_cols(g)
    nrows = len(Cidx)
    rp = [rank_modp(cols, nrows, p) for p in P]
    s = f"{name}: rank mod p {P[0]} = {rp[0]}, mod {P[1]} = {rp[1]}"
    re = None
    if EXACT:
        re = rank_exact(cols, nrows); s += f", EXACT over Z[i] = {re}"
    print(s + f"   ({time.time()-t0:.1f}s)", flush=True)
    return rp, re

print(f"m={m} f={f}: dim W = {dimW}, dim N = {dimN}, dim notN = {len(Cidx)}", flush=True)
D = {z: dir2(z) for z in range(1, n)}
res = {}
res['dir e1'] = fval('f(dir e1)', D[1])
if ALLB:
    for z in range(2, n): fval(f'f(dir {z})', D[z])
B = matmul(D[1], D[2])
res['B2'] = fval('f(dir e1 dir e2)', B)
if ALLB and m == 3:
    for z1, z2 in itertools.combinations(range(1, n), 2):
        if (z1, z2) == (1, 2): continue
        fval(f'f(dir {z1} dir {z2})', matmul(D[z1], D[z2]))
K = D[1]
for j in range(1, m): K = matmul(K, D[1 << j])
res['K'] = fval('f(kernel)', K)
# H = kernel of the wedge map V (x) Lambda^k V -> Lambda^(k+1) V has dimension n C(n,k) - C(n,k+1) and contains N,
# so f(anything) <= dim H - dim N:
from math import comb
print(f"dim H - dim N = {n*comb(n,k) - comb(n,k+1) - dimN}  (upper bound for every f)")
c1 = res['dir e1'][1] if EXACT else None
c2 = res['B2'][1] if EXACT else None
dk = res['K'][1] if EXACT else min(res['K'][0])
print(f"SUMMARY m={m} f={f}: cap={c1}, f(B2)={c2}, dem={dk}", flush=True)
