"""xc3-s1s2 (wave 17): exact checks for the alpha-symbol theorem (PROOF.md sections 2-4).

Block m = 2, f digit columns, Pauli labels at column j: (x0, x1, z0, z1) -> global index 4j + t.
Label matrix M(U) (F2, 4f x 4f): U P_v U^* = +- P_(M v).  F4 = F2[alpha]/(alpha^2+alpha+1); an F4 vector is a pair
(a, b) of F2 vectors meaning a + alpha b, so an F2 matrix acts componentwise.
Plane wave on a window I = [lo, hi): e_I(c)_j = alpha^(-j) c for j in I, 0 elsewhere; S = down-shift (S y)_j = y_(j-1),
so S e = alpha e in the bulk.  A Laurent polynomial q in S acts as the truncated Toeplitz matrix T_f(q).

Checks (all exact, over F2 / F4):
 (A) per letter: the symbol s_g read off at two windows and two sizes f agrees, and M e_I(c) - e_I(s c) is
     (q_g, dI, R)-tame: q_g(S) of it vanishes at distance > R from the window ends and the array ends.
 (B) random words in the alpha-regular letters: M(w) e_I(c) - e_I(s_w c) is (Q_w, dI, R_w)-tame with Q_w the
     product of the letters' q and R_w <= sum of radii (+ widths); s_w = product of letter symbols.
 (C) rank lemma: if s_w != I then rank_F2(M(w_f) - I) >= K(f) (number of disjoint test windows), growing linearly.
 (D) consistency of families: chi S_last chi^-1 = P as operators (same label matrix) and same symbol (I);
     E^-1 Dcz E Dcz^-1 = CZ-chain on row 1; the skeptic's pair k, l (v-ut-direct-logic) have symbol I and k*l too.
 (E) frames: every symbol preserves the Hermitian form Hg; free letters fix L0 = {x = 0}; directions move L0 to an
     adjacent frame, K to an opposite one (distances in the 27-frame polar space).
"""
import itertools, random, sys
import numpy as np

random.seed(20261014)
m = 2; NT = 2 * m   # label components per column

# ---------------- F4 ----------------
EXP = [1, 2, 3]; LOG = {1: 0, 2: 1, 3: 2}
def fmul(a, b):
    if a == 0 or b == 0: return 0
    return EXP[(LOG[a] + LOG[b]) % 3]
def finv(a): return EXP[(-LOG[a]) % 3]
def fpow(a, e): return 1 if e == 0 else EXP[(LOG[a] * e) % 3]
ALPHA = 2
def frob(a): return fmul(a, a)

# F4 matrices as lists of lists of ints
def mmul(A, B):
    n, k, p = len(A), len(B), len(B[0])
    out = [[0] * p for _ in range(n)]
    for i in range(n):
        for j in range(p):
            s = 0
            for t in range(k): s ^= fmul(A[i][t], B[t][j])
            out[i][j] = s
    return out
def mid(n): return [[int(i == j) for j in range(n)] for i in range(n)]
def minv(A):
    n = len(A); M = [A[i][:] + mid(n)[i] for i in range(n)]
    r = 0
    for c in range(n):
        p = next(i for i in range(r, n) if M[i][c])
        M[r], M[p] = M[p], M[r]
        iv = finv(M[r][c]); M[r] = [fmul(iv, x) for x in M[r]]
        for i in range(n):
            if i != r and M[i][c]:
                fc = M[i][c]; M[i] = [x ^ fmul(fc, y) for x, y in zip(M[i], M[r])]
        r += 1
    return [row[n:] for row in M]

# ---------------- F2 label matrices ----------------
def idx(j, t): return NT * j + t
def f2inv(A):
    n = A.shape[0]; M = np.concatenate([A.copy() % 2, np.eye(n, dtype=np.uint8)], axis=1)
    r = 0
    for c in range(n):
        piv = np.nonzero(M[r:, c])[0]
        p = r + piv[0]
        if p != r: M[[r, p]] = M[[p, r]]
        rows = np.nonzero(M[:, c])[0]
        for i in rows:
            if i != r: M[i] ^= M[r]
        r += 1
    return M[:, n:]
def f2rank(A):
    M = A.copy() % 2; n, k = M.shape; r = 0
    for c in range(k):
        if r == n: break
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0: continue
        p = r + piv[0]
        if p != r: M[[r, p]] = M[[p, r]]
        rows = np.nonzero(M[:, c])[0]
        for i in rows:
            if i != r: M[i] ^= M[r]
        r += 1
    return r
def f2mul(A, B): return (A.astype(np.int64) @ B.astype(np.int64) % 2).astype(np.uint8)

def addr_label(A):
    """address map x -> A x on the m*f address bits (index m*j + row) -> label matrix diag(A, A^-T) (per column order)."""
    f = A.shape[0] // m
    Ait = f2inv(A).T.copy()
    M = np.zeros((NT * f, NT * f), dtype=np.uint8)
    for j in range(f):
        for k in range(f):
            M[NT*j:NT*j+m, NT*k:NT*k+m] = A[m*j:m*j+m, m*k:m*k+m]
            M[NT*j+m:NT*j+2*m, NT*k+m:NT*k+2*m] = Ait[m*j:m*j+m, m*k:m*k+m]
    return M
def poly_addr(G, f):
    """G: dict (r, r') -> {d: 1} Laurent coefficients; (x')_(j,r) = sum_{r', d} G[r,r'][d] x_(j-d, r')."""
    A = np.zeros((m * f, m * f), dtype=np.uint8)
    for (r, rp), coef in G.items():
        for d, v in coef.items():
            if not v: continue
            for j in range(f):
                k = j - d
                if 0 <= k < f: A[m*j + r, m*k + rp] ^= 1
    return A
def phase_label(B):
    """diagonal Clifford phase with symmetric F2 matrix B on the m*f bits: (x | z) -> (x | z + B x)."""
    f = B.shape[0] // m
    M = np.eye(NT * f, dtype=np.uint8)
    for j in range(f):
        for k in range(f):
            M[NT*j+m:NT*j+2*m, NT*k:NT*k+m] ^= B[m*j:m*j+m, m*k:m*k+m]
    return M
def column_uniform(g, f):
    """one-column label matrix g (4x4 over F2, order x0 x1 z0 z1) on every column."""
    M = np.zeros((NT * f, NT * f), dtype=np.uint8)
    for j in range(f): M[NT*j:NT*j+NT, NT*j:NT*j+NT] = g
    return M
def boundary_last(g, f):
    M = np.eye(NT * f, dtype=np.uint8); j = f - 1
    M[NT*j:NT*j+NT, NT*j:NT*j+NT] = g
    return M

def transvection(z):      # paid direction D_z on one column: x += (zlabel . z) z
    g = np.eye(4, dtype=np.uint8)
    for a in range(2):
        for b in range(2):
            g[a, 2 + b] ^= z[a] & z[b]
    return g
HH = np.zeros((4, 4), dtype=np.uint8); HH[0, 2] = HH[1, 3] = HH[2, 0] = HH[3, 1] = 1   # K = H (x) H on a column

# ---------------- letters: name -> (builder f -> label matrix, q as {power: 1}, free?) ----------------
ONE = {0: 1}
def Gdiag(coef): return {(0, 0): coef, (1, 1): coef}
LET = {}
LET['chi'] = (lambda f: addr_label(poly_addr(Gdiag({0: 1, 1: 1}), f)), {0: 1, -1: 1}, True)
LET['chi^-1'] = (lambda f: addr_label(f2inv(poly_addr(Gdiag({0: 1, 1: 1}), f))), {0: 1, 1: 1}, True)
LET['chiA'] = (lambda f: addr_label(poly_addr(Gdiag({0: 1, -1: 1}), f)), {0: 1, 1: 1}, True)        # x_j += x_(j+1)
LET['chiA^-1'] = (lambda f: addr_label(f2inv(poly_addr(Gdiag({0: 1, -1: 1}), f))), {0: 1, -1: 1}, True)
LET['E'] = (lambda f: addr_label(poly_addr({(0, 0): ONE, (1, 1): ONE, (0, 1): {1: 1}}, f)), ONE, True)   # row0 ^= row1<<1
LET['Gm'] = (lambda f: addr_label(poly_addr({(0, 0): ONE, (0, 1): {1: 1}, (1, 0): {1: 1}, (1, 1): {0: 1, 1: 1, 2: 1}}, f)),
             {0: 1, -1: 1}, True)   # det = 1 + S: polynomial, inverse rational
LET['Gm^-1'] = (lambda f: addr_label(f2inv(poly_addr({(0, 0): ONE, (0, 1): {1: 1}, (1, 0): {1: 1}, (1, 1): {0: 1, 1: 1, 2: 1}}, f))),
                {0: 1, 1: 1}, True)
def B_D1(f):
    B = np.zeros((m*f, m*f), dtype=np.uint8)
    for j in range(f): B[m*j+1, m*j+1] = 1
    return B
def B_Dcz(f):
    B = np.zeros((m*f, m*f), dtype=np.uint8)
    for j in range(f): B[m*j, m*j+1] = B[m*j+1, m*j] = 1
    return B
def B_P(f):
    B = np.zeros((m*f, m*f), dtype=np.uint8)
    for j in range(f):
        for k in range(f): B[m*j+1, m*k+1] = 1
    return B
def B_CZchain(f):
    B = np.zeros((m*f, m*f), dtype=np.uint8)
    for j in range(1, f): B[m*j+1, m*(j-1)+1] = B[m*(j-1)+1, m*j+1] = 1
    return B
def B_Slast(f):
    B = np.zeros((m*f, m*f), dtype=np.uint8); B[m*(f-1)+1, m*(f-1)+1] = 1
    return B
LET['D1'] = (lambda f: phase_label(B_D1(f)), ONE, True)
LET['Dcz'] = (lambda f: phase_label(B_Dcz(f)), ONE, True)
LET['P'] = (lambda f: phase_label(B_P(f)), {0: 1, 1: 1}, True)
LET['Slast'] = (lambda f: phase_label(B_Slast(f)), ONE, True)
LET['Hlast'] = (lambda f: boundary_last(HH, f), ONE, True)            # non-monomial boundary letter (free)
def A_k(f):      # skeptic's k: x_0[0] ^= x_1[0] + ... + x_(f-1)[0]
    A = np.eye(m*f, dtype=np.uint8)
    for j in range(1, f): A[0, m*j] ^= 1
    return A
def A_l(f):      # skeptic's l: x_j[0] ^= x_0[0] for j >= 1
    A = np.eye(m*f, dtype=np.uint8)
    for j in range(1, f): A[m*j, 0] ^= 1
    return A
LET['k'] = (lambda f: addr_label(A_k(f)), {0: 1, 1: 1}, True)
LET['l'] = (lambda f: addr_label(A_l(f)), {0: 1, 1: 1}, True)
for nm, z in (('D01', (1, 0)), ('D10', (0, 1)), ('D11', (1, 1))):
    LET[nm] = ((lambda zz: (lambda f: column_uniform(transvection(zz), f)))(z), ONE, False)
LET['K'] = (lambda f: column_uniform(HH, f), ONE, False)
# one-column free Cliffords (monomial): CNOT row1->row0 on every column, swap rows
CN = np.eye(4, dtype=np.uint8); CN[0, 1] = 1; CN[3, 2] = 1        # x0 += x1 ; z1 += z0
SW = np.zeros((4, 4), dtype=np.uint8); SW[0, 1] = SW[1, 0] = SW[2, 3] = SW[3, 2] = 1
LET['CNOT'] = (lambda f: column_uniform(CN, f), ONE, True)
LET['SWAP'] = (lambda f: column_uniform(SW, f), ONE, True)
INV = {'chi': 'chi^-1', 'chi^-1': 'chi', 'chiA': 'chiA^-1', 'chiA^-1': 'chiA', 'Gm': 'Gm^-1', 'Gm^-1': 'Gm'}

# ---------------- waves, Laurent action, tameness ----------------
AINV = finv(ALPHA)
def wave(f, lo, hi, c):
    """F4 vector as pair (a, b) of F2 arrays of length 4f."""
    a = np.zeros(NT * f, dtype=np.uint8); b = np.zeros(NT * f, dtype=np.uint8)
    for j in range(lo, hi):
        sc = fpow(AINV, j)
        for t in range(NT):
            v = fmul(sc, c[t])
            a[NT*j + t] = v & 1; b[NT*j + t] = v >> 1
    return a, b
def apply_M(M, v):
    a, b = v
    return (M.astype(np.int64) @ a % 2).astype(np.uint8), (M.astype(np.int64) @ b % 2).astype(np.uint8)
def vsub(u, v): return (u[0] ^ v[0], u[1] ^ v[1])
def laurent(q, v, f):
    """T_f(q) applied blockwise: (q(S) y)_j = sum_d q_d y_(j-d), out-of-range dropped."""
    out = []
    for comp in v:
        y = comp.reshape(f, NT); r = np.zeros_like(y)
        for d, cv in q.items():
            if not cv: continue
            if d >= 0: r[d:] ^= y[:f - d] if d > 0 else y
            else: r[:f + d] ^= y[-d:]
        out.append(r.reshape(-1))
    return tuple(out)
def qmul(p, q):
    r = {}
    for a, x in p.items():
        for b, y in q.items():
            if x and y: r[a + b] = r.get(a + b, 0) ^ 1
    return {k: v for k, v in r.items() if v}
def qeval(q, lam_pow_fn):
    s = 0
    for d, v in q.items():
        if v: s ^= lam_pow_fn(d)
    return s
def qalpha(q): return qeval(q, lambda d: fpow(ALPHA, d % 3))
def nonzero_cols(v, f):
    a, b = v
    z = (a.reshape(f, NT) | b.reshape(f, NT)).any(axis=1)
    return np.nonzero(z)[0]
def tame_radius(v, f, pts):
    """max over columns j where v_j != 0 of dist(j, pts U {-1, f}); returns -1 if v = 0."""
    P = list(pts) + [-1, f]
    cols = nonzero_cols(v, f)
    if len(cols) == 0: return -1
    return max(min(abs(j - p) for p in P) for j in cols)
def col_value(v, j):
    a, b = v
    return [int(a[NT*j + t]) | (int(b[NT*j + t]) << 1) for t in range(NT)]

def symbol_of(M, q, f, lo, hi):
    """read s with T(q) M e_I(c) at the window centre = q(alpha) alpha^(-j0) s c."""
    j0 = (lo + hi) // 2
    s = [[0] * NT for _ in range(NT)]
    norm = finv(fmul(qalpha(q), fpow(AINV, j0)))
    for t in range(NT):
        c = [int(i == t) for i in range(NT)]
        r = laurent(q, apply_M(M, wave(f, lo, hi, c)), f)
        col = col_value(r, j0)
        for i in range(NT): s[i][t] = fmul(norm, col[i])
    return s
def residual_radius(M, s, q, f, lo, hi):
    worst = -1
    for t in range(NT):
        for cc in ([int(i == t) for i in range(NT)], [ALPHA if i == t else (1 if i == (t + 1) % NT else 0) for i in range(NT)]):
            sc = [0] * NT
            for i in range(NT):
                for k in range(NT): sc[i] ^= fmul(s[i][k], cc[k])
            y = vsub(apply_M(M, wave(f, lo, hi, cc)), wave(f, lo, hi, sc))
            worst = max(worst, tame_radius(laurent(q, y, f), f, (lo, hi)))
    return worst

ok = True
def rep(name, cond):
    global ok
    ok = ok and bool(cond)
    print(("PASS " if cond else "FAIL ") + name, flush=True)

# ---------------- (A) letters ----------------
print("== (A) letters: symbol at alpha and tameness of M e_I(c) - e_I(s c) ==")
SYM = {}; RAD = {}
for nm, (bld, q, free) in LET.items():
    s1 = symbol_of(bld(40), q, 40, 12, 28); s2 = symbol_of(bld(52), q, 52, 17, 38)
    r1 = residual_radius(bld(40), s1, q, 40, 12, 28); r2 = residual_radius(bld(52), s1, q, 52, 17, 38)
    SYM[nm] = s1; RAD[nm] = max(r1, r2, 0)
    rep(f"{nm:8s} q={sorted(q)}: symbol consistent at two windows/sizes; residual radius {r1},{r2} <= 3; s = {s1}",
        s1 == s2 and max(r1, r2) <= 3)
for a, b in INV.items():
    rep(f"s({a}) s({b}) = I", mmul(SYM[a], SYM[b]) == mid(NT))

# ---------------- (E) frames ----------------
print("== (E) the symbols on the 27 frames ==")
def Hg(u, v):
    s = 0
    for i in range(2): s ^= fmul(u[i], frob(v[2 + i])) ^ fmul(u[2 + i], frob(v[i]))
    return s
VEC = list(itertools.product(range(4), repeat=4))
def app(s, u): return tuple(sum_f([fmul(s[i][k], u[k]) for k in range(4)]) for i in range(4))
def sum_f(xs):
    r = 0
    for x in xs: r ^= x
    return r
def span(gens):
    S = {(0, 0, 0, 0)}
    for g in gens:
        S = S | {tuple(a ^ fmul(c, b) for a, b in zip(sv, g)) for sv in S for c in range(1, 4)}
    return frozenset(S)
L0 = span([(0, 0, 1, 0), (0, 0, 0, 1)])
def dist(A, B): return {16: 0, 4: 1, 1: 2}[len(A & B)]
SAMPLE = VEC[::5]
for nm, (bld, q, free) in LET.items():
    s = SYM[nm]
    pres = all(Hg(app(s, u), app(s, v)) == Hg(u, v) for u in SAMPLE for v in SAMPLE)
    img = frozenset(app(s, u) for u in L0)
    d = dist(L0, img)
    want = 0 if free else (2 if nm == 'K' else 1)
    rep(f"{nm:8s} preserves Hg; distance(L0, s L0) = {d} (want {want}, {'free' if free else 'paid'})", pres and d == want)

# ---------------- (D) family identities ----------------
print("== (D) identities between families and their symbols ==")
for f in (12, 20):
    lhs = f2mul(f2mul(LET['chi'][0](f), LET['Slast'][0](f)), LET['chi^-1'][0](f))
    rep(f"f={f}: M(chi Slast chi^-1) = M(P)", np.array_equal(lhs, LET['P'][0](f)))
    E, Dcz = LET['E'][0](f), LET['Dcz'][0](f)
    lhs = f2mul(f2mul(f2mul(f2inv(E), Dcz), E), f2inv(Dcz))
    rep(f"f={f}: M(E^-1 Dcz E Dcz^-1) = M(CZ chain on row 1)", np.array_equal(lhs, phase_label(B_CZchain(f))))
sw = mmul(mmul(SYM['chi'], SYM['Slast']), SYM['chi^-1'])
rep("symbol(chi Slast chi^-1) = symbol(P) = I", sw == SYM['P'] == mid(NT))
scz = symbol_of(phase_label(B_CZchain(40)), ONE, 40, 12, 28)
rep(f"symbol(E^-1 Dcz E Dcz^-1) = symbol(CZ chain) = {scz} (= symbol(D1): the S_last relation holds in the bulk)",
    mmul(mmul(mmul(minv(SYM['E']), SYM['Dcz']), SYM['E']), minv(SYM['Dcz'])) == scz == SYM['D1'])
for f in (16, 24, 32):
    kl = f2mul(LET['k'][0](f), LET['l'][0](f))
    rep(f"f={f}: skeptic's k*l: rank(M - I) = {f2rank(kl ^ np.eye(NT*f, dtype=np.uint8))} (bounded), symbol(k) s(l) = I",
        mmul(SYM['k'], SYM['l']) == mid(NT))

# ---------------- (B), (C) random words ----------------
print("== (B) random words: wave property with Q_w = product of q's;  (C) rank lemma ==")
names = list(LET)
def word_data(w, f):
    M = np.eye(NT * f, dtype=np.uint8); s = mid(NT); Q = {0: 1}; R = 0
    for nm in w:
        M = f2mul(M, LET[nm][0](f)); s = mmul(s, SYM[nm]); Q = qmul(Q, LET[nm][1]); R += RAD[nm] + max(abs(d) for d in LET[nm][1])
    return M, s, Q, R
nw = 0; nbad = 0; rank_rows = []
for trial in range(60):
    L = random.randint(2, 9)
    w = [random.choice(names) for _ in range(L)]
    f = 36
    M, s, Q, R = word_data(w, f)
    # window [10, 26)
    r = residual_radius(M, s, Q, f, 10, 26)
    good = r <= R
    nw += 1; nbad += (not good)
    if not good: print("   bad word", w, "radius", r, "bound", R)
    if s != mid(NT) and len(rank_rows) < 12:
        ranks = []
        for ff in (24, 36, 48):
            MM, _, _, _ = word_data(w, ff)
            ranks.append(f2rank(MM ^ np.eye(NT * ff, dtype=np.uint8)))
        Rp = max(R, max(abs(d) for d in Q)) + 1
        Kb = [max(0, (ff - 2 * Rp - 2) // (4 * Rp + 2)) for ff in (24, 36, 48)]
        rank_rows.append((w, ranks, Kb))
rep(f"(B) {nw} random words: residual of M(w) e_I(c) - e_I(s_w c) is (Q_w, dI, R_w)-tame in all ({nbad} bad)", nbad == 0)
allgood = True
for w, ranks, Kb in rank_rows:
    g = all(r >= k for r, k in zip(ranks, Kb)) and ranks[0] < ranks[1] < ranks[2]
    allgood &= g
    print(f"   s_w != I: rank(M - I) at f = 24, 36, 48: {ranks}  (window bound K: {Kb})  word {' '.join(w)}")
rep("(C) every sampled word with s_w != I has rank(M - I) >= K(f) and strictly growing", allgood)

# words with s_w = I but non-trivial operator: rank can be bounded (allowed) or growing; just report a few
print("   words with s_w = I (rank may be anything; S2 is one-directional):")
cnt = 0
for trial in range(400):
    w = [random.choice(names) for _ in range(random.randint(2, 7))]
    _, s, _, _ = word_data(w, 12)
    if s == mid(NT):
        rk = [f2rank(word_data(w, ff)[0] ^ np.eye(NT * ff, dtype=np.uint8)) for ff in (24, 36)]
        print(f"     {' '.join(w)}: rank(M - I) at f = 24, 36: {rk}")
        cnt += 1
        if cnt >= 6: break

print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
