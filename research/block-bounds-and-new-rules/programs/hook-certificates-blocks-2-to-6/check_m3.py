# pf-logging, 13 Oct 2026.  INDEPENDENT CHECKER for the block-3 foundation (14W/5): shares no code with gen_m3.py.
# Everything is rebuilt from the definitions:
#   V = C^8 with basis e_x, x in F_2^3 (integers 0..7, xor);  X_z e_x = e_{x xor z};
#   dir z = (1/2)((1+i) Id + (1-i) X_z)  (the engine's direction letter);  K = dir 1 . dir 2 . dir 4;
#   W = V (x) Lambda^4 V, g acts as g (x) Lambda^4 g, (Lambda^4 g)[U,T] = det g[U,T] (U, T sorted 4-sets);
#   S = {(x,T): x in T},  N = span{e_s : s in S}  (dim 280),  f(g) := dim(N + gN) - dim N = rank g_W[S^c, S].
# Checks:
#   (a) cap: for every z != 0, 200 explicit vectors u in N with (2 dir z)_W u in N, exactly in Z[i], with pairwise
#       disjoint supports  ->  rank (dir z)_W[S^c,S] <= 280 - 200 = 80.
#   (b) dem: cert_m3_dem.json names a square submatrix of K_W[S^c,S]; its determinant mod p is nonzero
#       -> rank K_W[S^c,S] over Q(i) >= its size (a nonzero image mod p forces a nonzero determinant).
#   (c) B2 = dir 1 . dir 2: cert_m3_B2.json lists exact Z[i] vectors v in N; check (2dir1 2dir2)_W v in N exactly and
#       that the vectors are independent (rank mod p = their number) -> f(B2) <= 280 - #vectors.
#   (d) M_0 (Lemma C, alpha = 0 block): cert_m3_M0.json minor nonsingular mod p -> rank M_0 >= its size.
# Usage: python3 -I check_m3.py   (reads the three JSON files from the current folder)
import itertools, json, sys

n, k = 8, 4
popc = lambda x: bin(x).count('1')
SUBS = list(itertools.combinations(range(n), k))
inS = lambda x, T: x in T
NS = sum(1 for x in range(n) for T in SUBS if inS(x, T))
assert NS == 280
ok_all = True
def report(name, ok, msg):
    global ok_all
    ok_all &= ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)

# ---------- exact Gaussian integers (pairs) ----------
def mul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def add(a, b): return (a[0] + b[0], a[1] + b[1])
ZERO, ONE = (0, 0), (1, 0)
def matmul(A, B):
    return [[(lambda r, c: __import__('functools').reduce(add, [mul(A[r][t], B[t][c]) for t in range(n)], ZERO))(r, c)
             for c in range(n)] for r in range(n)]
def two_dir(z):                    # 2 dir z = (1+i) Id + (1-i) X_z, exact
    M = [[ZERO] * n for _ in range(n)]
    for x in range(n):
        M[x][x] = add(M[x][x], (1, 1))
        M[x ^ z][x] = add(M[x ^ z][x], (1, -1))
    return M
def det_exact(A):                  # Laplace expansion along the first row (different routine from the generator)
    if len(A) == 1: return A[0][0]
    tot = ZERO
    for j in range(len(A)):
        if A[0][j] == ZERO: continue
        minor = [row[:j] + row[j + 1:] for row in A[1:]]
        t = mul(A[0][j], det_exact(minor))
        tot = add(tot, t if j % 2 == 0 else (-t[0], -t[1]))
    return tot
def wedge_sign_sorted(seq):        # e_{s1} ^ ... ^ e_{sk} = sign * e_{sorted}, or 0 if repeated
    if len(set(seq)) < len(seq): return 0, None
    inv = sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] > seq[j])
    return (-1 if inv % 2 else 1), tuple(sorted(seq))
_WEDGE_CACHE = {}
def apply_W(g, vec):               # g_W applied to a sparse vector {(x,T): gaussian int}; exact
    out = {}
    for (x, T), c in vec.items():
        if c == ZERO: continue
        ck = (id(g), T)            # (Lambda^4 g) e_T = sum_U det g[U,T] e_U  (cached per T)
        if ck not in _WEDGE_CACHE:
            col = {}
            for U in SUBS:
                d = det_exact([[g[u][t] for t in T] for u in U])
                if d != ZERO: col[U] = d
            _WEDGE_CACHE[ck] = (g, col)
        col = _WEDGE_CACHE[ck][1]
        for y in range(n):
            if g[y][x] == ZERO: continue
            gy = mul(g[y][x], c)
            for U, d in col.items():
                key = (y, U); out[key] = add(out.get(key, ZERO), mul(gy, d))
    return {kk: v for kk, v in out.items() if v != ZERO}

# ---------- (a) cap for every direction ----------
for z in range(1, n):
    g = two_dir(z)
    U_vecs = []
    for a in range(n):
        b = a ^ z
        if a > b: continue
        for w in itertools.combinations(range(n), k - 1):
            u = {}
            for x, sg in ((a, 1), (b, -1)):
                s, srt = wedge_sign_sorted((x,) + w)       # e_x (x) (e_x ^ e_w)
                if s: u[(x, srt)] = (sg * s, 0)
            if u: U_vecs.append(u)
    # exactness: (2 dir z)_W u has support in S
    bad = sum(1 for u in U_vecs if any(not inS(y, Uu) for (y, Uu) in apply_W(g, u)))
    supp = [set(u) for u in U_vecs]
    allkeys = [key for s_ in supp for key in s_]
    disjoint = len(allkeys) == len(set(allkeys)) and all(inS(x, T) for (x, T) in allkeys)
    report(f"cap dir {z}", bad == 0 and disjoint and len(U_vecs) == 200,
           f"{len(U_vecs)} vectors in N, images outside N: {bad}, disjoint supports: {disjoint} -> f(dir {z}) <= {NS - len(U_vecs)}")

# ---------- mod p tools ----------
def is_prime(q):
    if q < 2: return False
    d = 2
    while d * d <= q:
        if q % d == 0: return False
        d += 1
    return True
def det_modp(A, p):
    A = [r[:] for r in A]; N = len(A); det = 1
    for c in range(N):
        pr = next((r for r in range(c, N) if A[r][c] % p), None)
        if pr is None: return 0
        if pr != c: A[c], A[pr] = A[pr], A[c]; det = -det
        det = det * A[c][c] % p
        iv = pow(A[c][c], p - 2, p)
        for r in range(c + 1, N):
            if A[r][c] % p:
                f = A[r][c] * iv % p
                A[r] = [(A[r][j] - f * A[c][j]) % p for j in range(N)]
    return det % p
def rank_modp(A, p):
    A = [r[:] for r in A]; R = len(A); C = len(A[0]); rk = 0
    for c in range(C):
        pr = next((r for r in range(rk, R) if A[r][c] % p), None)
        if pr is None: continue
        A[rk], A[pr] = A[pr], A[rk]
        iv = pow(A[rk][c], p - 2, p)
        for r in range(rk + 1, R):
            if A[r][c] % p:
                f = A[r][c] * iv % p
                A[r] = [(A[r][j] - f * A[rk][j]) % p for j in range(C)]
        rk += 1
    return rk

def K_modp(p, I):
    h = pow(2, p - 2, p)
    def dirp(z):
        M = [[0] * n for _ in range(n)]
        for x in range(n):
            M[x][x] = (M[x][x] + (1 + I) * h) % p
            M[x ^ z][x] = (M[x ^ z][x] + (1 - I) * h) % p
        return M
    def mm(A, B): return [[sum(A[r][t] * B[t][c] for t in range(n)) % p for c in range(n)] for r in range(n)]
    return mm(mm(dirp(1), dirp(2)), dirp(4))

# ---------- (b) dem ----------
C = json.load(open("cert_m3_dem.json"))
p, I = C["p"], C["I"]
rows = [(r[0], tuple(r[1])) for r in C["rows"]]; cols = [(c[0], tuple(c[1])) for c in C["cols"]]
okfmt = (is_prime(p) and (I * I + 1) % p == 0 and len(rows) == len(cols) == len(set(rows)) == len(set(cols))
         and all(len(U) == k and list(U) == sorted(set(U)) and not inS(y, U) for y, U in rows)
         and all(len(T) == k and list(T) == sorted(set(T)) and inS(x, T) for x, T in cols))
K = K_modp(p, I)
def Kw(y, U, x, T):
    return K[y][x] * det_modp([[K[u][t] for t in T] for u in U], p) % p
Mn = [[Kw(y, U, x, T) for (x, T) in cols] for (y, U) in rows]
d = det_modp(Mn, p)
report("dem m=3 f=5", okfmt and d != 0, f"minor of K_W[S^c,S] of size {len(rows)}, det mod {p} = {d} -> dem >= {len(rows)}")
dem3 = len(rows) if (okfmt and d) else 0

# ---------- (c) B2 ----------
C2 = json.load(open("cert_m3_B2.json"))
basis = [(b[0], tuple(b[1])) for b in C2["basis_S"]]
okb = len(basis) == len(set(basis)) == NS and all(inS(x, T) for x, T in basis)
g2 = matmul(two_dir(1), two_dir(2))
vecs = [{basis[int(key)]: (c[0], c[1]) for key, c in v.items()} for v in C2["vectors"]]
bad = 0
for v in vecs:
    if any(not inS(y, U) for (y, U) in apply_W(g2, v)): bad += 1
pB, IB = 1000000007 + 2, None                     # independent prime choice for the independence test
while not (is_prime(pB) and pB % 4 == 1): pB += 4 if pB % 4 == 1 else 1
IB = next(pow(gg, (pB - 1) // 4, pB) for gg in range(2, 500) if pow(gg, (pB - 1) // 2, pB) == pB - 1)
bidx = {b: i for i, b in enumerate(basis)}
Vm = [[0] * NS for _ in vecs]
for r, v in enumerate(vecs):
    for key, c in v.items(): Vm[r][bidx[key]] = (c[0] + c[1] * IB) % pB
rk = rank_modp(Vm, pB)
report("f(dir1 dir2) m=3 f=5", okb and bad == 0 and rk == len(vecs),
       f"{len(vecs)} exact vectors v in N with B2_W v in N (failures {bad}), independent (rank mod {pB} = {rk}) -> f(B2) <= {NS - rk}")

# ---------- (d) M_0 block ----------
C3 = json.load(open("cert_m3_M0.json"))
p3, I3 = C3["p"], C3["I"]
r3 = [tuple(U) for U in C3["rows"]]; c3 = [tuple(T) for T in C3["cols"]]
okf3 = (is_prime(p3) and (I3 * I3 + 1) % p3 == 0 and len(r3) == len(c3) == len(set(r3)) == len(set(c3))
        and all(0 not in U and len(U) == k and list(U) == sorted(set(U)) for U in r3)
        and all(0 in T and len(T) == k and list(T) == sorted(set(T)) for T in c3))
K3 = K_modp(p3, I3)
def M0(U, T):                      # sum_s [K_W X_s]_{(0,U),(0,T)} = sum_s K[0,s] det K[U, (t xor s)_t]
    return sum(K3[0][s] * det_modp([[K3[u][t ^ s] for t in T] for u in U], p3) for s in range(n)) % p3
d3 = det_modp([[M0(U, T) for T in c3] for U in r3], p3)
report("M_0 block m=3 f=5", okf3 and d3 != 0, f"minor of size {len(r3)}, det mod {p3} = {d3} -> rank M_0 >= {len(r3)}; 8 x = {8*len(r3)}")

from fractions import Fraction
print(f"SUMMARY m=3 f=5: cap <= {NS-200}, dem >= {dem3}, f(B2) <= {NS-rk}: bound {Fraction(dem3, NS-200)} W; all checks {'PASS' if ok_all else 'FAIL'}")
sys.exit(0 if ok_all else 1)
