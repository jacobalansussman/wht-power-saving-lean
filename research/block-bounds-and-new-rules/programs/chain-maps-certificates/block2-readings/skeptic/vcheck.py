"""v-xc3-s1s2-compute (wave 17 skeptic): independent recomputation of xc3-s1s2's Theorem S / Theorem B2 inputs.

Own code, written from the definitions in xc3-s1s2/PROOF.md sections 1-2 (not copied from s1s2_check.py).
Conventions: bits b[j][r] (column j = 0..f-1, row r = 0,1); label vector per column block (x_j0, x_j1, z_j0, z_j1),
global index 4j + t.  Address map U|b> = |A b>  ->  X^a -> X^(A a), Z^c -> Z^(A^-T c).  Phase z += B x (B symmetric).
S is the down-shift (S y)_j = y_(j-1); T_f(p) truncated; wave e_I(c)_j = lam^-j c.
F4 = {0,1,2=alpha,3=alpha^2}; an F4 vector is stored as a pair of F2 vectors (v = v0 + alpha v1).

Parts:
 (A) regularity by definition: find q in a candidate list with T_f(q) M deep block rows banded Toeplitz (radius <= 6),
     same n_d for every deep row and for f in {23, 31, 40}; symbol s = n(alpha)/q(alpha).  Controls singular at alpha
     (R3 = row ^= row<<3, G3b) must FAIL.
 (B) Hg preserved by every symbol; 27 Hg-frames; free letters fix L0; directions move L0 one step; K two steps.
 (C) words: s(w) = product; direct read of s(w) from M(w) e_I(c) at f = 160 through T(Q_w) at the window centre;
     rank_F2(M(w_f) - I) at f = 16, 24, 32, 40, 48; every word with s(w) != I must have rank >= K(f) and growing.
 (D) the 27 Fermat lines: reduction mod 2 = all 27 generators of Hstd, delta over Q(omega) (mod 7 and mod 13) equals
     the generator distance on all 729 pairs; an explicit isometry Hg -> Hstd found by search, mapping frames to frames.
"""
import itertools, random, sys
import numpy as np

random.seed(20261014)
# ---------------- F4 ----------------
MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
INV = [None, 1, 3, 2]
ALPHA = 2
def fpow(a, e):
    r = 1
    e %= 3
    for _ in range(e): r = MUL[r][a]
    return r
def frob(a): return MUL[a][a]

def m4mul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    C = [[0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = 0
            for t in range(k): s ^= MUL[A[i][t]][B[t][j]]
            C[i][j] = s
    return C
def m4id(n=4): return [[int(i == j) for j in range(n)] for i in range(n)]
def m4inv(A):
    n = len(A); M = [row[:] + m4id(n)[i] for i, row in enumerate(A)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c]), None)
        if p is None: return None
        M[c], M[p] = M[p], M[c]
        iv = INV[M[c][c]]; M[c] = [MUL[iv][x] for x in M[c]]
        for r in range(n):
            if r != c and M[r][c]:
                fct = M[r][c]; M[r] = [x ^ MUL[fct][y] for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]
def m4vec(A, v):
    out = []
    for i in range(len(A)):
        x = 0
        for t in range(len(v)): x ^= MUL[A[i][t]][v[t]]
        out.append(x)
    return out

# ---------------- F2 matrices ----------------
def f2mul(A, B): return (A.astype(np.float64) @ B.astype(np.float64) % 2).astype(np.uint8)
def f2rank(A):
    A = A.copy() % 2; r = 0; rows, cols = A.shape
    for c in range(cols):
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0: continue
        p = r + piv[0]
        if p != r: A[[r, p]] = A[[p, r]]
        nz = np.nonzero(A[:, c])[0]
        nz = nz[nz != r]
        A[nz] ^= A[r]
        r += 1
        if r == rows: break
    return r
def f2inv(A):
    n = A.shape[0]; M = np.concatenate([A % 2, np.eye(n, dtype=np.uint8)], 1)
    for c in range(n):
        piv = np.nonzero(M[c:, c])[0]
        assert len(piv), "singular"
        p = c + piv[0]
        if p != c: M[[c, p]] = M[[p, c]]
        nz = np.nonzero(M[:, c])[0]; nz = nz[nz != c]
        M[nz] ^= M[c]
    return M[:, n:].copy()

def bidx(j, r): return 2 * j + r
def lidx_x(j, r): return 4 * j + r
def lidx_z(j, r): return 4 * j + 2 + r

def label_from_addr(A, f):
    """A: 2f x 2f bit map. Returns the 4f x 4f label matrix in per-column ordering."""
    AinvT = f2inv(A).T
    M = np.zeros((4 * f, 4 * f), np.uint8)
    for j in range(f):
        for r in range(2):
            for jj in range(f):
                for rr in range(2):
                    M[lidx_x(jj, rr), lidx_x(j, r)] = A[bidx(jj, rr), bidx(j, r)]
                    M[lidx_z(jj, rr), lidx_z(j, r)] = AinvT[bidx(jj, rr), bidx(j, r)]
    return M
def label_from_phase(B, f):
    M = np.eye(4 * f, dtype=np.uint8)
    for j in range(f):
        for r in range(2):
            for jj in range(f):
                for rr in range(2):
                    if B[bidx(jj, rr), bidx(j, r)]:
                        M[lidx_z(jj, rr), lidx_x(j, r)] ^= 1
    return M
def column_uniform(g, f):
    M = np.zeros((4 * f, 4 * f), np.uint8)
    for j in range(f): M[4 * j:4 * j + 4, 4 * j:4 * j + 4] = g
    return M
def on_column(g, f, col):
    M = np.eye(4 * f, dtype=np.uint8); M[4 * col:4 * col + 4, 4 * col:4 * col + 4] = g; return M

def addr_poly(G, f):
    """G: dict (r_out, r_in) -> list of shifts d (coefficient 1 at S^d); b'_(j,ro) = sum b_(j-d, ri)."""
    A = np.zeros((2 * f, 2 * f), np.uint8)
    for (ro, ri), ds in G.items():
        for d in ds:
            for j in range(f):
                if 0 <= j - d < f: A[bidx(j, ro), bidx(j - d, ri)] ^= 1
    return A

def symplectic_ok(M, f):
    J = np.zeros((4 * f, 4 * f), np.uint8)
    for j in range(f):
        for r in range(2):
            J[lidx_x(j, r), lidx_z(j, r)] = 1; J[lidx_z(j, r), lidx_x(j, r)] = 1
    return np.array_equal(f2mul(f2mul(M.T, J), M), J)

# one-column constants (x0, x1, z0, z1)
def g_from_perm_swap(pairs):
    g = np.eye(4, dtype=np.uint8)
    for a, b in pairs:
        g[[a, b]] = g[[b, a]]
    return g
G_SWAPQ = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]], np.uint8)       # swap the two qubits
G_CNOT = np.array([[1,0,0,0],[1,1,0,0],[0,0,1,1],[0,0,0,1]], np.uint8)        # b1 ^= b0: x1 += x0; z0 += z1
G_H0 = g_from_perm_swap([(0, 2)])
G_H1 = g_from_perm_swap([(1, 3)])
G_HH = g_from_perm_swap([(0, 2), (1, 3)])
def g_transvection(v):      # x += (z . v) v   (xc-round2 / s1s2 convention for a paid direction)
    g = np.eye(4, dtype=np.uint8)
    for i in range(2):
        for k in range(2):
            if v[i] and v[k]: g[i, 2 + k] ^= 1
    return g
G_KX = np.array([[1,0,1,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]], np.uint8)          # x += z
G_H11 = (lambda C: f2mul(f2mul(C, G_H0), C))(G_CNOT)                            # Hadamard along (1,1)-ish
G_SS = np.array([[1,0,0,0],[0,1,0,0],[1,0,1,0],[0,1,0,1]], np.uint8)          # S on both qubits (z += x)

# ---------------- letters ----------------
def full_B(f, pairs_fn):
    B = np.zeros((2 * f, 2 * f), np.uint8)
    for (a, b) in pairs_fn(f):
        B[a, b] ^= 1
        if a != b: B[b, a] ^= 1
    return B

LETTERS = {}
def L(name, kind, fn): LETTERS[name] = (kind, fn)
L("chi", "free", lambda f: label_from_addr(addr_poly({(0, 0): [0, 1], (1, 1): [0, 1]}, f), f))
L("chi_ac", "free", lambda f: label_from_addr(addr_poly({(0, 0): [0, -1], (1, 1): [0, -1]}, f), f))
L("chi_row1", "free", lambda f: label_from_addr(addr_poly({(0, 0): [0], (1, 1): [0, 1]}, f), f))
L("E", "free", lambda f: label_from_addr(addr_poly({(0, 0): [0], (0, 1): [1], (1, 1): [0]}, f), f))
L("Gm", "free", lambda f: label_from_addr(addr_poly({(0, 0): [0], (0, 1): [1], (1, 0): [1], (1, 1): [0, 1, 2]}, f), f))
L("G5", "free", lambda f: label_from_addr(addr_poly({(0, 0): [0, 2], (1, 1): [0, 1, 2, 3], (0, 1): [1]}, f), f))  # G = [[1+S^2, S],[0,(1+S)^3]], det (1+S)^5
L("D1", "free", lambda f: label_from_phase(full_B(f, lambda f: [(bidx(j, r), bidx(j, r)) for j in range(f) for r in range(2)]), f))
L("Dcz", "free", lambda f: label_from_phase(full_B(f, lambda f: [(bidx(j, 0), bidx(j, 1)) for j in range(f)]), f))
L("CZch", "free", lambda f: label_from_phase(full_B(f, lambda f: [(bidx(j, 0), bidx(j + 1, 1)) for j in range(f - 1)]), f))
L("P", "free", lambda f: label_from_phase(full_B(f, lambda f: [(bidx(j, 1), bidx(k, 1)) for j in range(f) for k in range(j, f)]), f))
L("Slast", "free", lambda f: label_from_phase(full_B(f, lambda f: [(bidx(f - 1, 0), bidx(f - 1, 0))]), f))
L("Hlast", "free", lambda f: on_column(G_H0, f, f - 1))
L("Hfirst", "free", lambda f: on_column(G_H1, f, 0))
def A_k(f):
    A = np.eye(2 * f, dtype=np.uint8)
    for j in range(1, f): A[bidx(0, 0), bidx(j, 0)] = 1
    return A
def A_l(f):
    A = np.eye(2 * f, dtype=np.uint8)
    for j in range(1, f): A[bidx(j, 0), bidx(0, 0)] = 1
    return A
L("k", "free", lambda f: label_from_addr(A_k(f), f))
L("l", "free", lambda f: label_from_addr(A_l(f), f))
L("CNOT", "free", lambda f: column_uniform(G_CNOT, f))
L("SWAPQ", "free", lambda f: column_uniform(G_SWAPQ, f))
L("SS", "free", lambda f: column_uniform(G_SS, f))
L("Dt01", "paid1", lambda f: column_uniform(g_transvection((0, 1)), f))
L("Dt10", "paid1", lambda f: column_uniform(g_transvection((1, 0)), f))
L("Dt11", "paid1", lambda f: column_uniform(g_transvection((1, 1)), f))
L("H0", "paid1", lambda f: column_uniform(G_H0, f))
L("H1", "paid1", lambda f: column_uniform(G_H1, f))
L("H11", "paid1", lambda f: column_uniform(G_H11, f))
L("HH", "target", lambda f: column_uniform(G_HH, f))
L("KX", "target", lambda f: column_uniform(G_KX, f))
# controls singular at alpha
L("R3", "control", lambda f: label_from_addr(addr_poly({(0, 0): [0, 3], (1, 1): [0, 3]}, f), f))
L("G3b", "control", lambda f: label_from_addr(addr_poly({(0, 0): [0, 1], (0, 1): [1], (1, 0): [1], (1, 1): [0]}, f), f))

# ---------------- Laurent operators ----------------
def T_scalar(q, f):
    """q: dict d -> 1 (F2). Returns T_f(q) (x) I4 as 4f x 4f."""
    T = np.zeros((4 * f, 4 * f), np.uint8)
    for d, c in q.items():
        if not c: continue
        for j in range(f):
            if 0 <= j - d < f:
                for t in range(4): T[4 * j + t, 4 * (j - d) + t] ^= 1
    return T
def qmul(p, q):
    r = {}
    for a, x in p.items():
        for b, y in q.items():
            if x and y: r[a + b] = r.get(a + b, 0) ^ 1
    return {d: 1 for d, c in r.items() if c}
def qeval(q, lam_pow):
    s = 0
    for d, c in q.items():
        if c: s ^= lam_pow(d)
    return s
def qwidth(q): return max(abs(d) for d in q) if q else 0

Q_CANDS = []
for tot in range(0, 11):
    for a in range(0, 6):
        b = tot - a
        if not (0 <= b <= 5): continue
        q = {0: 1}
        for _ in range(a): q = qmul(q, {0: 1, 1: 1})
        for _ in range(b): q = qmul(q, {0: 1, -1: 1})
        Q_CANDS.append((f"(1+S)^{a}(1+S^-1)^{b}", q))

def deep_toeplitz(M, q, f, r=10):
    """If deep block rows (dist to {-1,f} > r) of T_f(q) M are banded (|d| <= r) Toeplitz with one n, return n."""
    QM = f2mul(T_scalar(q, f), M)
    n = None
    for j in range(f):
        if not (j + 1 > r and f - j > r): continue
        row = QM[4 * j:4 * j + 4, :]
        nz_cols = np.nonzero(row.any(0))[0]
        if len(nz_cols) and (abs(nz_cols // 4 - j)).max() > r: return None
        nd = {}
        for d in range(-r, r + 1):
            k = j - d
            if 0 <= k < f:
                blk = row[:, 4 * k:4 * k + 4]
                if blk.any(): nd[d] = blk.copy()
        if n is None: n = nd
        else:
            if set(n) != set(nd) or any(not np.array_equal(n[d], nd[d]) for d in n): return None
    return n

def symbol_at(n, q, lam_pow):
    S = [[0] * 4 for _ in range(4)]
    for d, blk in n.items():
        lp = lam_pow(d)
        for i in range(4):
            for k in range(4):
                if blk[i, k]: S[i][k] ^= lp
    qv = qeval(q, lam_pow)
    if qv == 0: return None
    iq = INV[qv]
    return [[MUL[iq][x] for x in row] for row in S]

alpha_pow = lambda d: fpow(ALPHA, d)
FS = [23, 31, 40]
ok_all = True
def rep(name, cond, extra=""):
    global ok_all
    if not cond and not name.startswith("CONTROL"): ok_all = False
    print(("PASS " if cond else "FAIL ") + name + (("  " + extra) if extra else ""))
    sys.stdout.flush()

print("== (A) regularity by definition and symbols at alpha ==")
SYM = {}; QOF = {}; RAD = {}
for name, (kind, fn) in LETTERS.items():
    found = None
    for qn, q in Q_CANDS:
        if qeval(q, alpha_pow) == 0: continue
        ns = []
        for f in FS:
            M = fn(f)
            assert symplectic_ok(M, f), (name, f)
            n = deep_toeplitz(M, q, f)
            if n is None: break
            ns.append(n)
        if len(ns) == len(FS):
            same = all(set(ns[0]) == set(x) and all(np.array_equal(ns[0][d], x[d]) for d in ns[0]) for x in ns)
            if same:
                found = (qn, q, ns[0]); break
    if kind == "control":
        rep(f"CONTROL {name} singular at alpha: no regular data found (expected)", found is None)
        continue
    if found is None:
        rep(f"{name}: regular at alpha", False); continue
    qn, q, n = found
    s = symbol_at(n, q, alpha_pow)
    SYM[name] = s; QOF[name] = q; RAD[name] = max([qwidth(q)] + [abs(d) for d in n])
    rep(f"{name} ({kind}): regular, q = {qn}, band {RAD[name]}", s is not None and m4inv(s) is not None,
        "symbol rows " + str(s))

# inverse letters: own regular data for M^-1
for name in list(SYM):
    kind, fn = LETTERS[name]
    found = None
    for qn, q in Q_CANDS:
        if qeval(q, alpha_pow) == 0: continue
        ns = []
        for f in FS:
            n = deep_toeplitz(f2inv(fn(f)), q, f)
            if n is None: break
            ns.append(n)
        if len(ns) == len(FS) and all(set(ns[0]) == set(x) and all(np.array_equal(ns[0][d], x[d]) for d in ns[0]) for x in ns):
            found = (qn, q, ns[0]); break
    if found is None:
        rep(f"{name}^-1 regular", False); continue
    qn, q, n = found
    si = symbol_at(n, q, alpha_pow)
    SYM[name + "^-1"] = si; QOF[name + "^-1"] = q; RAD[name + "^-1"] = max([qwidth(q)] + [abs(d) for d in n])
    rep(f"{name}^-1 regular (q = {qn}) and s(k^-1) = s(k)^-1", m4mul(SYM[name], si) == m4id())

print("\n== (B) Hermitian form, frames, moves ==")
def Hg(u, v):
    s = 0
    for i in range(2): s ^= MUL[u[i]][frob(v[2 + i])] ^ MUL[u[2 + i]][frob(v[i])]
    return s
VEC = list(itertools.product(range(4), repeat=4))
def app(s, u): return tuple(int(x) for x in m4vec(s, list(u)))
for name, s in SYM.items():
    BAS = [tuple(int(i == t) for i in range(4)) for t in range(4)]
    good = all(Hg(app(s, u), app(s, v)) == Hg(u, v) for u in BAS for v in BAS)
    good = good and all(Hg(app(s, u), app(s, v)) == Hg(u, v) for u in VEC[::5] for v in VEC[::3])
    rep(f"{name}: symbol preserves Hg", good)
def span(gens):
    S = {(0, 0, 0, 0)}
    for g in gens:
        S = S | {tuple(a ^ MUL[c][b] for a, b in zip(x, g)) for x in S for c in range(1, 4)}
    return frozenset(S)
FR = set()
NZ = [v for v in VEC if any(v)]
for u in NZ:
    if Hg(u, u): continue
    for v in NZ:
        if Hg(v, v) or Hg(u, v): continue
        S2 = span([u, v])
        if len(S2) == 16: FR.add(S2)
FR = sorted(FR, key=lambda S: sorted(S))
rep(f"27 Hg-frames (got {len(FR)})", len(FR) == 27)
L0 = span([(0, 0, 1, 0), (0, 0, 0, 1)])
rep("L0 = {x = 0} is a frame", L0 in FR)
def dist(A, B): return {16: 0, 4: 1, 1: 2}[len(A & B)]
def img(s, S): return frozenset(app(s, u) for u in S)
for name, s in SYM.items():
    base = name.replace("^-1", "")
    kind = LETTERS[base][0]
    d = dist(L0, img(s, L0))
    want = {"free": 0, "paid1": 1, "target": 2}[kind]
    rep(f"{name} ({kind}): dist(L0, s L0) = {d}", d == want and img(s, L0) in FR)

print("\n== (C) words: product symbol, direct read, rank growth ==")
NAMES = [n for n in SYM]
CACHE = {}
def letter_M(nm, f):
    key = (nm, f)
    if key not in CACHE:
        base = nm.replace("^-1", "")
        X = LETTERS[base][1](f)
        if nm.endswith("^-1"): X = f2inv(X)
        CACHE[key] = X
    return CACHE[key]
def word_M(w, f):
    M = np.eye(4 * f, dtype=np.uint8)
    for nm in w: M = f2mul(M, letter_M(nm, f))
    return M
def word_sym(w):
    s = m4id()
    for nm in w: s = m4mul(s, SYM[nm])
    return s
def wave(f, lo, hi, c):    # c: F4 4-vector -> pair of F2 vectors
    v0 = np.zeros(4 * f, np.uint8); v1 = np.zeros(4 * f, np.uint8)
    for j in range(lo, hi):
        lp = fpow(ALPHA, -j)
        for t in range(4):
            e = MUL[lp][c[t]]
            v0[4 * j + t] = e & 1; v1[4 * j + t] = e >> 1
    return v0, v1
def direct_symbol(w, f=160, lo=40, hi=121):
    M = word_M(w, f)
    Q = {0: 1}
    for nm in w: Q = qmul(Q, QOF[nm])
    TQ = T_scalar(Q, f)
    TM = f2mul(TQ, M)
    jc = (lo + hi) // 2
    cols = []
    for t in range(4):
        c = [int(i == t) for i in range(4)]
        v0, v1 = wave(f, lo, hi, c)
        y0 = TM.astype(np.int64) @ v0 % 2; y1 = TM.astype(np.int64) @ v1 % 2
        blk = [int(y0[4 * jc + i]) | (int(y1[4 * jc + i]) << 1) for i in range(4)]
        # blk = Q(alpha) alpha^-jc s c  ->  s c = blk * alpha^jc / Q(alpha)
        qa = qeval(Q, alpha_pow)
        if qa == 0: return None
        fac = MUL[fpow(ALPHA, jc)][INV[qa]]
        cols.append([MUL[fac][x] for x in blk])
    return [[cols[k][i] for k in range(4)] for i in range(4)]
WORDS = []
free_names = [n for n in NAMES if LETTERS[n.replace("^-1", "")][0] == "free"]
for _ in range(70):
    t = random.randint(2, 8)
    WORDS.append([random.choice(NAMES) for _ in range(t)])
# targeted words
WORDS += [["chi", "Slast", "chi^-1"], ["chi", "Slast", "chi^-1", "P^-1"], ["k", "l"], ["E^-1", "Dcz", "E", "Dcz^-1"],
          ["chi", "Hlast", "chi^-1"], ["Gm", "Dt01", "Gm^-1", "Dt01^-1"], ["chi", "Dt11", "chi^-1", "Dt11^-1"],
          ["E", "H0", "E^-1", "H0^-1"], ["P", "Dt10", "P^-1", "Dt10^-1"], ["l", "H0", "l^-1", "H0^-1"],
          ["k", "Dt01", "k^-1", "Dt01^-1"], ["chi_ac", "chi"], ["Gm", "E^-1", "Gm^-1", "E"]]
FW = [16, 24, 32, 40, 48]
n_bad_direct = 0; n_bad_rank = 0; summary = []
for w in WORDS:
    sw = word_sym(w)
    sd = direct_symbol(w)
    if sd != sw: n_bad_direct += 1; print("   direct read mismatch:", w, sd, sw)
    ranks = [f2rank(word_M(w, f) ^ np.eye(4 * f, dtype=np.uint8)) for f in FW]
    R = sum(2 * RAD[nm] for nm in w); Rp = R + 1
    Kf = [max(0, (f - 4 * Rp - 1) // (4 * Rp + 2) + 1) for f in FW]
    isI = sw == m4id()
    if not isI:
        if any(r < k for r, k in zip(ranks, Kf)) or not (ranks[-1] > ranks[0]): n_bad_rank += 1; print("   RANK FAIL:", w, ranks, Kf)
    summary.append((isI, ranks, w))
rep(f"direct read of s(w) from M(w) equals the product of letter symbols on {len(WORDS)} words", n_bad_direct == 0)
rep(f"every word with s(w) != I has rank(M - I) >= K(f) and growing ({sum(1 for s in summary if not s[0])} words)", n_bad_rank == 0)
bounded = [s for s in summary if s[0]]
print(f"   words with s(w) = I: {len(bounded)}; their rank sequences:")
for isI, ranks, w in bounded: print("    ", ranks, " ".join(w))
print("   sample of s(w) != I:")
for isI, ranks, w in [s for s in summary if not s[0]][:8]: print("    ", ranks, " ".join(w))

print("\n== (D) the 27 Fermat lines: exact over Z[omega], reduction mod 2, distances ==")
# Z[omega] as pairs (a, b) = a + b w, w^2 = -1 - w
def zadd(x, y): return (x[0] + y[0], x[1] + y[1])
def zneg(x): return (-x[0], -x[1])
def zmul(x, y):
    a, b = x; c, d = y
    # (a + b w)(c + d w) = ac + (ad + bc) w + bd w^2 = (ac - bd) + (ad + bc - bd) w
    return (a * c - b * d, a * d + b * c - b * d)
ZERO, ONE, W = (0, 0), (1, 0), (0, 1)
def wpow(i):
    r = ONE
    for _ in range(i % 3): r = zmul(r, W)
    return r
def zdet(M):
    n = len(M)
    if n == 1: return M[0][0]
    s = ZERO
    for c in range(n):
        minor = [row[:c] + row[c + 1:] for row in M[1:]]
        t = zmul(M[0][c], zdet(minor))
        s = zadd(s, t if c % 2 == 0 else zneg(t))
    return s
def zrank(rows):
    n = len(rows); m = len(rows[0])
    for r in range(min(n, m), 0, -1):
        for R in itertools.combinations(range(n), r):
            for C in itertools.combinations(range(m), r):
                if zdet([[rows[i][j] for j in C] for i in R]) != ZERO: return r
    return 0
LINES = []
for (a, b), (c, d) in [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]:
    for i in range(3):
        for j in range(3):
            v1 = [ZERO] * 4; v1[b] = ONE; v1[a] = zneg(wpow(i))
            v2 = [ZERO] * 4; v2[d] = ONE; v2[c] = zneg(wpow(j))
            LINES.append((v1, v2))
def on_cubic(v):
    s = ZERO
    for x in v: s = zadd(s, zmul(x, zmul(x, x)))
    return s == ZERO
rep("27 lines lie on the Fermat cubic (generators and a few combinations)",
    all(on_cubic(v1) and on_cubic(v2) and on_cubic([zadd(x, zmul((2, 1), y)) for x, y in zip(v1, v2)]) for v1, v2 in LINES))
# reduce mod 2: omega -> alpha; a + b w -> (a mod 2) + (b mod 2) alpha
def red(x): return (x[0] % 2) ^ MUL[x[1] % 2][ALPHA]
def Hstd(u, v):
    s = 0
    for i in range(4): s ^= MUL[u[i]][frob(v[i])]
    return s
RL = [span([tuple(red(x) for x in v1), tuple(red(x) for x in v2)]) for v1, v2 in LINES]
rep("27 reductions are distinct 2-dim Hstd-isotropic subspaces", len(set(RL)) == 27 and all(len(S) == 16 for S in RL)
    and all(Hstd(u, v) == 0 for S in RL for u in S for v in S))
mism = 0; counts = {}
for U in range(27):
    for V in range(27):
        rk = zrank([list(LINES[U][0]), list(LINES[U][1]), list(LINES[V][0]), list(LINES[V][1])])
        delta = rk - 2
        dd = dist(RL[U], RL[V])
        counts[dd] = counts.get(dd, 0) + 1
        if delta != dd: mism += 1
rep(f"exact over Z[omega]: delta(Y_U, Y_V) = generator distance on all 729 ordered pairs (mismatches {mism}); by distance {counts}", mism == 0)
# isometry Phi: (F4^4, Hg) -> (F4^4, Hstd): columns p0, p1 (x-basis), q0, q1 (z-basis)
iso = [v for v in NZ if Hstd(v, v) == 0]
Phi = None
for p0 in iso:
    for q0 in iso:
        if Hstd(p0, q0) != 1: continue
        perp = [v for v in iso if Hstd(v, p0) == 0 and Hstd(v, q0) == 0]
        for p1 in perp:
            for q1 in perp:
                if Hstd(p1, q1) == 1:
                    Phi = [p0, p1, q0, q1]; break
            if Phi: break
        if Phi: break
    if Phi: break
def phi_app(u):
    out = [0, 0, 0, 0]
    for k in range(4):
        for i in range(4): out[i] ^= MUL[u[k]][Phi[k][i]]
    return tuple(out)
BAS = [tuple(int(i == t) for i in range(4)) for t in range(4)]
rep("Phi is an isometry Hg -> Hstd", Phi is not None and all(Hstd(phi_app(u), phi_app(v)) == Hg(u, v) for u in BAS for v in BAS))
img_frames = [frozenset(phi_app(u) for u in S) for S in FR]
rep("Phi maps the 27 Hg-frames onto the 27 reduced Fermat lines", set(img_frames) == set(RL))
print("\nALL CHECKS PASS" if ok_all else "\nSOME CHECK FAILED")
