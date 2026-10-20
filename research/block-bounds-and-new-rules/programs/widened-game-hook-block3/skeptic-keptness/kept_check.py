# v28-kept (wave 16 skeptic).  Independent check of the d-tof-cert hook certificate, written from scratch.
# Works in the engine's own data space at f = 5 digit columns: V^(x)5, V = C^8, a letter g acts as g^(x)5
# (Lean: matAct f M = digitProd (fun _ => M) *v x).  W = V (x) Lambda^4 V is realised INSIDE V^(x)5 as the vectors
# antisymmetric in slots 1..4; N = span{ e_x (x) Alt(e_T) : x in T }.  No code shared with d-tof-cert.
#   A. letters typed from Lean (BinaryFrames.lean: dir z = ca*1 + cb*shift z, kernel = littleC^(x)3,
#      Letters.lean: permMat, phaseMat), checks K = dir e1 dir e2 dir e3;
#   B. keptness of N, EXACT, by: all 7 shifts, all 168 additive bijections, all 16 phaseMat z c, the S_8 generators
#      (0 1) and (0 1 ... 7), Toffoli, CCZ, CS, T-phase-type integer diagonals; arbitrary diagonals by the
#      "one content multiset per basis vector" argument (checked);
#   C. exact Z[i] check of the 200 explicit vectors u with (2R)^(x)5 u in N, mod-p independence (3 primes);
#   D. f(dir z) for all 7 z, f(dir z1 dir z2) for all 42 ordered independent pairs, f(K) -- mod-p ranks, 3 primes.
# Usage: python3 -I kept_check.py
import sys, itertools, time
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np

n, F = 8, 5
t0 = time.time()
def say(*a):
    print(*a, f"[{time.time()-t0:.1f}s]", flush=True)

# ---------- Gaussian-integer 8x8 matrices as (re, im) int64 pairs ----------
def gmat(re, im=None):
    re = np.array(re, dtype=np.int64)
    return (re, np.zeros_like(re) if im is None else np.array(im, dtype=np.int64))
def gmul(A, B):
    return (A[0] @ B[0] - A[1] @ B[1], A[0] @ B[1] + A[1] @ B[0])
def shift(z):           # Lean: shift z i j = 1 iff i + z = j
    M = np.zeros((n, n), np.int64)
    for i in range(n): M[i, i ^ z] = 1
    return gmat(M)
I8 = gmat(np.eye(n, dtype=np.int64))
def dir2(z):            # 2 * dir z = (1+i) 1 + (1-i) shift z
    S = shift(z)[0]
    return gmat(np.eye(n, dtype=np.int64) + S, np.eye(n, dtype=np.int64) - S)
def dot(a, b): return bin(a & b).count('1') & 1
def littleC2_kron3():   # 8 * kernel = (2 littleC)^(x)3, digit j of the address = bit j
    re = np.zeros((n, n), np.int64); im = np.zeros((n, n), np.int64)
    for x in range(n):
        for y in range(n):
            # prod over bits of (1+i) if equal else (1-i)
            c = (1, 0)
            for j in range(3):
                a = (1, 1) if ((x >> j) & 1) == ((y >> j) & 1) else (1, -1)
                c = (c[0] * a[0] - c[1] * a[1], c[0] * a[1] + c[1] * a[0])
            re[x, y], im[x, y] = c
    return (re, im)
K8 = littleC2_kron3()
P = gmul(gmul(dir2(1), dir2(2)), dir2(4))
assert np.array_equal(P[0], K8[0]) and np.array_equal(P[1], K8[1])
say("A. 8*kernel (littleC^(x)3, Lean) == (2dir e1)(2dir e2)(2dir e3): True")
D = dir2(1); D2 = gmul(D, D); S1 = shift(1)
assert np.array_equal(D2[0], 4 * S1[0]) and not D2[1].any()
say("A. dir(e1)^2 == shift e1: True")

def perm_from_map(G):   # Lean permMat G x y = 1 iff G x = y
    M = np.zeros((n, n), np.int64)
    for x in range(n): M[x, G[x]] = 1
    return gmat(M)
def diag_int(d_re, d_im=None):
    return gmat(np.diag(d_re), None if d_im is None else np.diag(d_im))
def tint(e):            # i^e for e in Z/4 as (re, im)
    return [(1, 0), (0, 1), (-1, 0), (0, -1)][e % 4]
def phaseMat(z, c):     # diagonal i^[z.x + c], exponent in F_2 (xor)
    vals = [tint(dot(z, x) ^ c) for x in range(n)]
    return diag_int([v[0] for v in vals], [v[1] for v in vals])

# ---------- V^(x)5, apply g^(x)5 to a batch of vectors ----------
def terms_of(g):
    """split an 8x8 Gaussian-integer matrix into row-monomial terms: list of (col index per row, re coef, im coef)"""
    nzs = [[y for y in range(n) if g[0][x, y] or g[1][x, y]] for x in range(n)]
    out = []
    for k in range(max(len(r) for r in nzs)):
        col = np.array([r[k] if k < len(r) else 0 for r in nzs])
        cre = np.array([g[0][x, nzs[x][k]] if k < len(nzs[x]) else 0 for x in range(n)], np.int64)
        cim = np.array([g[1][x, nzs[x][k]] if k < len(nzs[x]) else 0 for x in range(n)], np.int64)
        out.append((col, cre, cim))
    return out
def apply_g5(g, Vre, Vim):
    """Vre, Vim: (k, 8,8,8,8,8) int64.  returns g^(x)5 applied, exactly: (g v)[.., x, ..] = sum_y g[x,y] v[.., y, ..]"""
    T = terms_of(g)
    for ax in range(1, F + 1):
        shp = [1] * (F + 1); shp[ax] = n
        Rr = np.zeros_like(Vre); Ri = np.zeros_like(Vim)
        for (col, cre, cim) in T:
            Ar = np.take(Vre, col, axis=ax); Ai = np.take(Vim, col, axis=ax)
            cr = cre.reshape(shp); ci = cim.reshape(shp)
            Rr += cr * Ar - ci * Ai; Ri += cr * Ai + ci * Ar
        Vre, Vim = Rr, Ri
    return Vre, Vim

S4 = list(itertools.combinations(range(n), 4))
def sign_of(seq):
    s = 1; seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]: s = -s
    return s
def alt_vec(x, T):      # e_x (x) sum_pi sgn(pi) e_{T pi}, T sorted
    v = np.zeros((n,) * F, np.int64)
    for pi in itertools.permutations(range(4)):
        v[(x,) + tuple(T[j] for j in pi)] = sign_of(pi)
    return v
Wbasis = [(x, T) for x in range(n) for T in S4]
Nbasis = [(x, T) for (x, T) in Wbasis if x in T]
Obasis = [(x, T) for (x, T) in Wbasis if x not in T]
say(f"dim W = {len(Wbasis)}, dim N = {len(Nbasis)}, dim W/N = {len(Obasis)}")
NV = np.stack([alt_vec(x, T) for (x, T) in Nbasis])          # (280, 8^5)
Z = np.zeros_like(NV)

# position masks
idx = np.indices((n,) * F).reshape(F, -1).T                   # all 5-tuples
inN_support = np.array([t[0] in t[1:] for t in idx]).reshape((n,) * F)
def antisym_ok(V):
    for (a, b) in [(1, 2), (2, 3), (3, 4)]:
        if not np.array_equal(np.swapaxes(V, a + 1, b + 1), -V):
            return False
    return True
def in_N(Vre, Vim):
    """exact: every vector antisymmetric in slots 1..4 and supported where y0 in {y1..y4}"""
    if (Vre[:, ~inN_support] != 0).any() or (Vim[:, ~inN_support] != 0).any(): return False
    return antisym_ok(Vre) and antisym_ok(Vim)
assert in_N(NV, Z)

# ---------- B. keptness ----------
lin = []
for cols in itertools.product(range(1, 8), repeat=3):
    img = [0] * 8
    for x in range(8):
        y = 0
        for k in range(3):
            if (x >> k) & 1: y ^= cols[k]
        img[x] = y
    if len(set(img)) == 8: lin.append(img)
assert len(lin) == 168
free = {}
for z in range(1, 8): free[f"shift {z}"] = shift(z)
for j, G in enumerate(lin): free[f"perm GL3 #{j}"] = perm_from_map(G)
for z in range(8):
    for c in (0, 1): free[f"phaseMat {z} {c}"] = phaseMat(z, c)
free["S8 gen (0 1)"] = perm_from_map([1, 0, 2, 3, 4, 5, 6, 7])
free["S8 gen 8-cycle"] = perm_from_map([(x + 1) % 8 for x in range(8)])
free["Toffoli x2^=x0x1"] = perm_from_map([x ^ (4 if (x & 3) == 3 else 0) for x in range(8)])
free["CCZ"] = diag_int([-1 if x == 7 else 1 for x in range(8)])
free["CS"] = diag_int([0 if (x & 3) == 3 else 1 for x in range(8)], [1 if (x & 3) == 3 else 0 for x in range(8)])
free["int diag (2,3,5,7,11,13,17,19)"] = diag_int([2, 3, 5, 7, 11, 13, 17, 19])
free["Z4 pattern random"] = (lambda q: diag_int([tint(e)[0] for e in q], [tint(e)[1] for e in q]))(
    np.random.default_rng(1).integers(0, 4, 8))
free["scalar i (phaseMat 0 1)"] = phaseMat(0, 1)
bad = []
for name, g in free.items():
    Rr, Ri = apply_g5(g, NV, Z)
    if not in_N(Rr, Ri): bad.append(name)
say(f"B. letters tested: {len(free)}; image of N not inside N for: {bad if bad else 'none'}")
# arbitrary diagonals: each N basis vector must be supported on tuples with ONE content multiset
one_multiset = all(len({tuple(sorted(t)) for t in zip(*np.nonzero(v))}) == 1 for v in NV)
say(f"B. every N basis vector supported on one multiset of addresses (=> eigenvector of EVERY diagonal^(x)5): {one_multiset}")

# ---------- mod-p tools ----------
def primes_1mod4(start, k):
    out = []; q = start
    while len(out) < k:
        q -= 1
        if q % 4 == 1 and all(q % d for d in range(2, int(q ** 0.5) + 1)): out.append(q)
    return out
PR = primes_1mod4(2 ** 27, 3)
def sqrtm1(p):
    for g in range(2, 200):
        if pow(g, (p - 1) // 2, p) == p - 1: return pow(g, (p - 1) // 4, p)
def rank_modp(M, p):
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
Opos = tuple(np.array([(x,) + T for (x, T) in Obasis]).T)
Npos = tuple(np.array([(x,) + T for (x, T) in Nbasis]).T)
def out_coords(Vre, Vim, p, I):
    """W-coordinates outside N (coefficient at the sorted position) reduced mod p, rows = vectors"""
    return (Vre[(slice(None),) + Opos] + I * (Vim[(slice(None),) + Opos] % p)) % p
def fval(gs, p, I):
    Vre, Vim = NV, Z
    for g in gs: Vre, Vim = apply_g5(g, Vre, Vim)
    assert antisym_ok(Vre) and antisym_ok(Vim)
    return rank_modp(out_coords(Vre % p, Vim % p, p, I), p)

# ---------- C. explicit vectors ----------
pairs = [(a, a ^ 1) for a in range(n) if a < (a ^ 1)]
U = []
for (a, b) in pairs:
    for w in itertools.combinations(range(n), 3):
        v = np.zeros((n,) * F, np.int64)
        for (xx, s) in ((a, 1), (b, -1)):
            if xx in w: continue
            # e_xx (x) Alt(e_xx (x) e_w1 (x) e_w2 (x) e_w3), built directly by antisymmetrising the 4-tuple
            tup = (xx,) + w
            for pi in itertools.permutations(range(4)):
                v[(xx,) + tuple(tup[j] for j in pi)] += s * sign_of(pi)
        if v.any(): U.append(v)
U = np.stack(U)
assert in_N(U, np.zeros_like(U))
Rr, Ri = apply_g5(dir2(1), U, np.zeros_like(U))
say(f"C. explicit vectors u: {len(U)}; all u in N: True; (2dir e1)^(x)5 u in N for all (exact Z[i]): {in_N(Rr, Ri)}")
for p in PR:
    I = sqrtm1(p)
    rU = rank_modp(np.stack([u.reshape(-1) for u in U]) % p, p)
    say(f"C. p={p}: rank of the u's = {rU}  -> cap(dir e1) <= {len(Nbasis) - rU} (exact, over Q(i))")

# ---------- D. f values mod p ----------
for p in PR[:2]:
    I = sqrtm1(p)
    fd = [fval([dir2(z)], p, I) for z in range(1, 8)]
    say(f"D. p={p}: f(dir z), z=1..7 (mod-p lower bounds): {fd}")
    fk = fval([dir2(1), dir2(2), dir2(4)], p, I)
    say(f"D. p={p}: f(K) (mod-p lower bound) = {fk}")
p = PR[0]; I = sqrtm1(p)
f2 = {}
for z1 in range(1, 8):
    for z2 in range(1, 8):
        if z1 == z2: continue
        f2[(z1, z2)] = fval([dir2(z2), dir2(z1)], p, I)   # block dir z1 * dir z2 applied: rightmost first
say(f"D. p={p}: f(dir z1 dir z2) over 42 ordered independent pairs: values {sorted(set(f2.values()))}")
fz = fval([dir2(1), dir2(1)], p, I)
say(f"D. p={p}: f(dir e1 dir e1) = f(shift) = {fz}")
