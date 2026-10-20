# pf-logging, 13 Oct 2026.  INDEPENDENT CHECKER for cert_r1.json (t514-tight's g2 table, engine game, block 3, Y = O_c).
# Shares no code with gen_r1.py.  From the definitions:
#   V = C^8, e_x (x in F_2^3); 2 dir z = (1+i) Id + (1-i) X_z;  a class representative (sigma, q) is the monomial matrix
#   b e_x = i^{q[x]} e_{sigma[x]};  g = (2 dir 1) b (2 dir 1);  W = V (x) Lambda^4 V, g_W[(y,U),(x,T)] = g[y,x] det g[U,T].
#   Y = O_c = coordinates (x,T) with x not in T and T not an affine plane (t1^t2^t3^t4 != 0); comp = the other 336.
# Checks, per class: every certificate vector v (supported on Y) has (2 dir 1)_W v and g_W v supported on Y, exactly in
#   Z[i]; the vectors are independent (rank mod a prime = their number).  And once: a square submatrix of
#   (2 dir 1)_W[comp, Y] with nonzero determinant mod p.
# Conclusion per class (exact_g2.py's definition g2 = rank[(2D)_W ; g_W][comp, Y] - rank (2D)_W[comp, Y]):
#   g2 <= (224 - #vectors) - (size of the minor).
# Usage: python3 -I check_r1.py cert_r1.json
import itertools, json, sys
n, k = 8, 4
C = json.load(open(sys.argv[1]))
def mul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def add(a, b): return (a[0] + b[0], a[1] + b[1])
Z = (0, 0)
def matmul(A, B):
    return [[__import__('functools').reduce(add, (mul(A[r][t], B[t][c]) for t in range(n)), Z) for c in range(n)] for r in range(n)]
def two_dir(z):
    M = [[Z] * n for _ in range(n)]
    for x in range(n): M[x][x] = add(M[x][x], (1, 1)); M[x ^ z][x] = add(M[x ^ z][x], (1, -1))
    return M
def monomial(sig, q):
    ipow = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    M = [[Z] * n for _ in range(n)]
    for x in range(n): M[sig[x]][x] = ipow[q[x] % 4]
    return M
def det(A):                        # Laplace along the first row
    if len(A) == 1: return A[0][0]
    tot = Z
    for j in range(len(A)):
        if A[0][j] == Z: continue
        t = mul(A[0][j], det([r[:j] + r[j + 1:] for r in A[1:]]))
        tot = add(tot, t if j % 2 == 0 else (-t[0], -t[1]))
    return tot
SUBS = list(itertools.combinations(range(n), k))
inY = lambda x, T: (x not in T) and (T[0] ^ T[1] ^ T[2] ^ T[3]) != 0
Y = [(x, T) for x in range(n) for T in SUBS if inY(x, T)]
comp = [(x, T) for x in range(n) for T in SUBS if not inY(x, T)]
okY = [(c[0], tuple(c[1])) for c in C["Y"]] == Y
print(f"|Y| = {len(Y)}, |comp| = {len(comp)}, certificate's Y list matches: {okY}", flush=True)
def compW(g):                      # (g_W)[comp, Y] as a dict per Y column: list over comp rows
    wedge = {T: {U: det([[g[u][t] for t in T] for u in U]) for U in SUBS} for T in set(T for _, T in Y)}
    return [[mul(g[y][x], wedge[T][U]) for (y, U) in comp] for (x, T) in Y]   # column-major: cols[j][i]
def images_in_Y(cols, v):         # is the comp-part of g_W v zero?  v: dict j -> gaussian int
    acc = [Z] * len(comp)
    for j, c in v.items():
        col = cols[j]
        for i in range(len(comp)):
            e = col[i]
            if e != Z: acc[i] = add(acc[i], mul(e, c))
    return all(a == Z for a in acc)
def is_prime(q):
    if q < 2: return False
    d = 2
    while d * d <= q:
        if q % d == 0: return False
        d += 1
    return True
def rank_modp(rows, p):
    A = [r[:] for r in rows]; rk = 0; ncol = len(A[0]) if A else 0
    for c in range(ncol):
        pr = next((r for r in range(rk, len(A)) if A[r][c] % p), None)
        if pr is None: continue
        A[rk], A[pr] = A[pr], A[rk]; iv = pow(A[rk][c], p - 2, p)
        for r in range(rk + 1, len(A)):
            if A[r][c] % p:
                f = A[r][c] * iv % p; A[r] = [(A[r][j] - f * A[rk][j]) % p for j in range(ncol)]
        rk += 1
    return rk
def det_modp(A, p):
    A = [r[:] for r in A]; N = len(A); d = 1
    for c in range(N):
        pr = next((r for r in range(c, N) if A[r][c] % p), None)
        if pr is None: return 0
        if pr != c: A[c], A[pr] = A[pr], A[c]; d = -d
        d = d * A[c][c] % p; iv = pow(A[c][c], p - 2, p)
        for r in range(c + 1, N):
            if A[r][c] % p:
                f = A[r][c] * iv % p; A[r] = [(A[r][j] - f * A[c][j]) % p for j in range(N)]
    return d % p
D2 = two_dir(1)
colsD = compW(D2)
# ---- the minor of (2D)_W[comp, Y] ----
mn = C["minor_BD"]; p, I = mn["p"], mn["I"]
rows = [(r[0], tuple(r[1])) for r in mn["rows"]]; cols = [(c[0], tuple(c[1])) for c in mn["cols"]]
cidx = {b: i for i, b in enumerate(comp)}; yidx = {b: j for j, b in enumerate(Y)}
okm = (is_prime(p) and (I * I + 1) % p == 0 and len(set(rows)) == len(rows) == len(set(cols)) == len(cols)
       and all(r in cidx for r in rows) and all(c in yidx for c in cols))
d = det_modp([[(colsD[yidx[c]][cidx[r]][0] + colsD[yidx[c]][cidx[r]][1] * I) % p for c in cols] for r in rows], p) if okm else 0
rBD = len(rows) if (okm and d) else 0
print(f"[{'PASS' if rBD else 'FAIL'}] minor of (2D)_W[comp,Y]: size {len(rows)}, det mod {p} = {d} -> rank >= {rBD}", flush=True)
allok = bool(rBD) and okY
q2 = 1000000009
I2 = next(pow(g_, (q2 - 1) // 4, q2) for g_ in range(2, 500) if pow(g_, (q2 - 1) // 2, q2) == q2 - 1)
table = []
for cl in C["classes"]:
    g = matmul(matmul(D2, monomial(cl["sigma"], cl["q"])), D2)
    colsG = compW(g)
    vecs = [{int(j): (c[0], c[1]) for j, c in v.items()} for v in cl["vectors"]]
    bad = sum(1 for v in vecs if not (images_in_Y(colsD, v) and images_in_Y(colsG, v)))
    rk = rank_modp([[((v.get(j, Z)[0] + v.get(j, Z)[1] * I2) % q2) for j in range(len(Y))] for v in vecs], q2)
    ok = bad == 0 and rk == len(vecs)
    g2 = len(Y) - rk - rBD
    allok &= ok; table.append(g2)
    print(f"[{'PASS' if ok else 'FAIL'}] class sigma={cl['sigma']} q={cl['q']}: {len(vecs)} vectors, exact failures {bad}, "
          f"independent {rk} -> g2 <= {g2}", flush=True)
print(f"RESULT: {'PASS' if allok else 'FAIL'}; certified upper bounds on g2 over the {len(table)} classes: {sorted(table)}")
sys.exit(0 if allok else 1)
