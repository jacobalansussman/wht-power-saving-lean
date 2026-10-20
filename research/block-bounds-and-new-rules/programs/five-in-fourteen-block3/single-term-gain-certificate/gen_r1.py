# pf-logging, 13 Oct 2026.  GENERATOR (untrusted) of certificates for t514-tight's exact g2 table (Lemma R1, engine
# game, block 3, Y = O_c), wave16/t514-tight/exact_g2.py, out_exact_g2_Oc.txt.
#   W = V (x) Lambda^4 V (m = 3, f = 5); Y = O_c = span{e_(x,T) : x not in T, T not an affine plane (xor of T != 0)}
#   D = dir(e1); for each of the 14 class representatives b = perm(sigma) diag(i^q):  g = (2D) b (2D);
#   g2(b) = rank[(2D)_W ; g_W] restricted to (comp rows, Y cols) - rank[(2D)_W restricted]  (exact_g2's definition)
# Upper bound certificate: exact vectors v in Y with (2D)_W v in Y and g_W v in Y; L of them, independent, give
#   rank[stacked] <= 224 - L, hence g2 <= 224 - L - rank[(2D)_W[comp, Y]] <= 224 - L - 80 once a nonsingular 80-minor of
#   (2D)_W[comp, Y] is exhibited (also emitted here).
# Method: null space mod p with i -> I and with i -> -I (RREF; the two give real and imaginary parts), rational
# reconstruction, exact re-verification in Z[i].  Pure Python.  Output: cert_r1.json
import itertools, json, sys
from fractions import Fraction as Fr
n, k = 8, 4
REPS = [([3, 2, 4, 5, 7, 6, 0, 1], [0, 0, 3, 1, 1, 3, 2, 2]), ([6, 7, 5, 4, 2, 3, 1, 0], [0, 1, 0, 3, 1, 0, 1, 2]),
        ([1, 5, 7, 3, 6, 2, 0, 4], [0, 0, 3, 1, 3, 1, 0, 0]), ([7, 4, 2, 1, 5, 6, 0, 3], [0, 0, 2, 0, 1, 3, 3, 3]),
        ([6, 0, 1, 7, 3, 5, 4, 2], [0, 0, 1, 1, 3, 1, 0, 2]), ([0, 4, 5, 1, 6, 2, 3, 7], [0, 3, 3, 0, 3, 0, 0, 3]),
        ([2, 1, 6, 5, 4, 7, 0, 3], [0, 3, 2, 1, 2, 3, 2, 3]), ([7, 0, 5, 2, 6, 1, 4, 3], [0, 3, 3, 0, 0, 3, 1, 2]),
        ([2, 6, 3, 7, 4, 0, 5, 1], [0, 1, 2, 3, 1, 0, 1, 0]), ([1, 4, 6, 3, 2, 7, 5, 0], [0, 2, 1, 3, 3, 1, 0, 2]),
        ([1, 2, 0, 3, 5, 6, 4, 7], [0, 3, 3, 2, 0, 3, 3, 2]), ([1, 2, 4, 7, 3, 0, 6, 5], [0, 0, 0, 2, 0, 0, 2, 0]),
        ([7, 0, 4, 3, 6, 1, 5, 2], [0, 2, 2, 0, 3, 3, 1, 1]), ([0, 2, 6, 4, 7, 5, 1, 3], [0, 0, 0, 2, 3, 3, 1, 3])]
P = 2147483629
I = next(pow(g, (P - 1) // 4, P) for g in range(2, 1000) if pow(g, (P - 1) // 2, P) == P - 1)
def gm(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def ga(a, b): return (a[0] + b[0], a[1] + b[1])
Z0 = (0, 0); IP = [(1, 0), (0, 1), (-1, 0), (0, -1)]
def mm(A, B):
    out = [[Z0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = Z0
            for l in range(n): s = ga(s, gm(A[i][l], B[l][j]))
            out[i][j] = s
    return out
def twoD(z):
    M = [[Z0] * n for _ in range(n)]
    for x in range(n): M[x][x] = ga(M[x][x], (1, 1)); M[x ^ z][x] = ga(M[x ^ z][x], (1, -1))
    return M
def mono(sig, q):
    M = [[Z0] * n for _ in range(n)]
    for x in range(n): M[sig[x]][x] = IP[q[x] % 4]
    return M
PERMS = []
for pm in itertools.permutations(range(k)):
    sg = 1
    for a in range(k):
        for b in range(a + 1, k):
            if pm[a] > pm[b]: sg = -sg
    PERMS.append((pm, sg))
def det4(M):
    tot = Z0
    for pm, sg in PERMS:
        t = (sg, 0)
        for r in range(k):
            t = gm(t, M[r][pm[r]])
            if t == Z0: break
        tot = ga(tot, t)
    return tot
SUBS = list(itertools.combinations(range(n), k))
basis = [(x, T) for x in range(n) for T in SUBS]
plane = lambda T: (T[0] ^ T[1] ^ T[2] ^ T[3]) == 0
Y = [b for b in basis if b[0] not in b[1] and not plane(b[1])]
comp = [b for b in basis if not (b[0] not in b[1] and not plane(b[1]))]
print(f"|Y| = {len(Y)}, |comp| = {len(comp)}", flush=True)
def W_block(g):            # g_W[comp, Y] as list of rows (exact)
    wedge = {}
    for T in set(T for (_, T) in Y):
        wedge[T] = {U: det4([[g[u][t] for t in T] for u in U]) for U in SUBS}
    return [[gm(g[y][x], wedge[T][U]) for (x, T) in Y] for (y, U) in comp]
def tomod(a, Iv): return (a[0] + a[1] * Iv) % P
import numpy as np
def nullspace_modp(rows, ncols):
    """RREF null-space basis mod P (numpy, int64; P < 2^31 so products fit)"""
    A = np.array(rows, dtype=np.int64) % P; piv = []; R = 0
    for c in range(ncols):
        nz = np.nonzero(A[R:, c])[0]
        if len(nz) == 0: continue
        pr = R + nz[0]
        if pr != R: A[[R, pr]] = A[[pr, R]]
        A[R] = A[R] * pow(int(A[R, c]), P - 2, P) % P
        col = A[:, c].copy(); col[R] = 0; nzr = np.nonzero(col)[0]
        if len(nzr): A[nzr] = (A[nzr] - (col[nzr, None] * A[R][None, :]) % P) % P
        piv.append(c); R += 1
        if R == A.shape[0]: break
    free = [c for c in range(ncols) if c not in piv]
    basis_ = []
    for fc in free:
        v = [0] * ncols; v[fc] = 1
        for r, pc in enumerate(piv): v[pc] = int((-A[r, fc]) % P)
        basis_.append(v)
    return basis_, free
def ratrec(a):
    a %= P
    r0, r1, s0, s1 = P, a, 0, 1
    while r1 * r1 * 2 > P:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or s1 * s1 * 2 > P: return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)
D2 = twoD(1)
BD = W_block(D2)
inv2 = pow(2, P - 2, P)
out = {"Y": [[x, list(T)] for x, T in Y], "classes": []}
# 80-minor of (2D)_W[comp, Y]
Mm = [[tomod(v, I) for v in row] for row in BD]
A = [r[:] for r in Mm]; pivs = []; used = [False] * len(A)
for c in range(len(Y)):
    pr = next((r for r in range(len(A)) if not used[r] and A[r][c]), None)
    if pr is None: continue
    used[pr] = True; pivs.append((pr, c)); iv = pow(A[pr][c], P - 2, P)
    for r in range(len(A)):
        if not used[r] and A[r][c]:
            f = A[r][c] * iv % P; A[r] = [(A[r][j] - f * A[pr][j]) % P for j in range(len(Y))]
print(f"rank mod p of (2D)_W[comp,Y] = {len(pivs)}", flush=True)
out["minor_BD"] = {"p": P, "I": I, "rows": [[comp[r][0], list(comp[r][1])] for r, c in pivs],
                   "cols": [[Y[c][0], list(Y[c][1])] for r, c in pivs]}
for ci, (sig, q) in enumerate(REPS):
    g = mm(mm(D2, mono(sig, q)), D2)
    Bg = W_block(g)
    stacked = BD + Bg
    nsA, freeA = nullspace_modp([[tomod(v, I) for v in row] for row in stacked], len(Y))
    nsB, freeB = nullspace_modp([[tomod(v, P - I) for v in row] for row in stacked], len(Y))
    assert freeA == freeB
    vecs = []
    for va, vb in zip(nsA, nsB):
        ent = []
        for xa, xb in zip(va, vb):
            re = (xa + xb) * inv2 % P; im = (xa - xb) * inv2 % P * pow(I, P - 2, P) % P
            ent.append((ratrec(re), ratrec(im)))
        assert all(e[0] is not None and e[1] is not None for e in ent), "reconstruction failed"
        from math import lcm
        L = 1
        for e in ent: L = lcm(L, e[0].denominator, e[1].denominator)
        iv_ = [(int(e[0] * L), int(e[1] * L)) for e in ent]
        # exact verification
        for row in stacked:
            s = Z0
            for a_, b_ in zip(row, iv_):
                if b_ != Z0 and a_ != Z0: s = ga(s, gm(a_, b_))
            assert s == Z0, "exact check failed"
        vecs.append({str(j): list(v) for j, v in enumerate(iv_) if v != Z0})
    g2 = len(Y) - len(vecs) - len(pivs)
    print(f"class {ci} sigma={sig} q={q}: {len(vecs)} exact vectors -> g2 <= {g2}", flush=True)
    out["classes"].append({"sigma": sig, "q": q, "vectors": vecs})
json.dump(out, open("cert_r1.json", "w"))
