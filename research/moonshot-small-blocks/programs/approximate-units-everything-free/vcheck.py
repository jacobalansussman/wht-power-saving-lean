"""v-xw-border/vcheck.py -- skeptic's own rebuild (does not import xw2lib/keptlib).

Block 2, f digit columns, address x = sum_j (r0_j + 2 r1_j) 4^j (same address convention as the game).
Dense 4^f x 4^f matrices built from scratch.  Checks:
 T1 R8 = pi- u2 - pi- u3 at exact Q, random fibrewise (Y (x) I, I (x) Y'), random full matrices.
 T2 Proposition S: R8-seeded family with random perturbations vs its exact twin at Q; boundedness at fibrewise Y;
    1/eps growth at full matrices.
 T3 commutant: a1, b1 (mirror chain on row 1), lamps; [D1, .] at fibrewise values; a1 order; [a1, b1];
    order of <L, L'> in GL_f(F2).
 T4 braid D_z S_z D_z = S_z D_z S_z at Q and at generic fibrewise values.
 T5 Fact A: tau(a1^n) = 2^(min(2^v, f) - f).
 T6 K = D1 D2 check, and D1 D2 = D2 D1 at fibrewise values.
"""
import numpy as np, itertools, sys

a, b = (1 + 1j) / 2, (1 - 1j) / 2
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 20261014)


def bits(f):
    n = 4 ** f
    x = np.arange(n)
    r0 = np.stack([(x >> (2 * j)) & 1 for j in range(f)], 1)
    r1 = np.stack([(x >> (2 * j + 1)) & 1 for j in range(f)], 1)
    return x, r0, r1


def addr(r0, r1):
    f = r0.shape[1]
    return sum((r0[:, j] << (2 * j)) | (r1[:, j] << (2 * j + 1)) for j in range(f))


def perm_mat(dest, phase=None):
    n = dest.size
    M = np.zeros((n, n), complex)
    ph = np.ones(n) if phase is None else phase
    M[dest, np.arange(n)] = ph          # e_x -> ph[x] e_dest[x]
    return M


def letters(f):
    x, r0, r1 = bits(f)
    L = {}
    L["X1"] = perm_mat(addr(1 - r0, r1))
    L["X2"] = perm_mat(addr(r0, 1 - r1))
    # chain map a1: row1 of col j ^= row1 of col j-1 (original values); mirror b1: col j ^= col j+1
    n1 = r1.copy(); n1[:, 1:] ^= r1[:, :-1]; L["a1"] = perm_mat(addr(r0, n1))
    m1 = r1.copy(); m1[:, :-1] ^= r1[:, 1:]; L["b1"] = perm_mat(addr(r0, m1))
    n0 = r0.copy(); n0[:, 1:] ^= r0[:, :-1]; L["a0"] = perm_mat(addr(n0, r1))
    fx = r1.copy(); fx[:, 0] ^= 1; L["fX2"] = perm_mat(addr(r0, fx))
    lx = r0.copy(); lx[:, f - 1] ^= 1; L["lX1"] = perm_mat(addr(lx, r1))
    L["lZ1"] = perm_mat(x, (-1.0) ** r0[:, f - 1])
    # S_z = i^(parity(z & col)) per column, product over columns
    for z in (1, 2, 3):
        par = ((r0 * (z & 1)) ^ (r1 * ((z >> 1) & 1))).sum(1)
        L["S%d" % z] = perm_mat(x, 1j ** par)
    # C01: row1 ^= row0 per column
    L["C01"] = perm_mat(addr(r0, r1 ^ r0))
    return L, x, r0, r1


def kron_pow(g, f):
    out = np.ones((1, 1), complex)
    for _ in range(f):
        out = np.kron(g, out)
    return out


def fib(f, Y, which):
    """embed an operator Y on the row-`which` bits (2^f dim, index sum_j bit_j 2^j) into the 4^f space"""
    x, r0, r1 = bits(f)
    i0 = (r0 << np.arange(f)).sum(1); i1 = (r1 << np.arange(f)).sum(1)
    n = 4 ** f
    M = np.zeros((n, n), complex)
    if which == 0:
        # M[x', x] = Y[i0(x'), i0(x)] if i1(x') == i1(x)
        same = i1[:, None] == i1[None, :]
        M = np.where(same, Y[i0[:, None], i0[None, :]], 0)
    else:
        same = i0[:, None] == i0[None, :]
        M = np.where(same, Y[i1[:, None], i1[None, :]], 0)
    return M


def exactD(f, z):
    """D_z = (a I + b X^z) on every column, built directly on the address"""
    x, r0, r1 = bits(f)
    n = 4 ** f
    M = np.eye(n, dtype=complex)
    for j in range(f):
        Xj = perm_mat(x ^ (z << (2 * j)))
        M = (a * np.eye(n) + b * Xj) @ M
    return M


def rnd(d):
    return (rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))) / np.sqrt(2 * d)


def rel(A, B):
    return np.linalg.norm(A - B) / max(np.linalg.norm(A), np.linalg.norm(B), 1e-300)


SQX = np.array([[a, b], [b, a]])
out = []
P = lambda *s: (print(*s), out.append(" ".join(str(t) for t in s)))

P("== T6: K = D1 D2 ; exact D1 = Q(x)I ; fibrewise commutation")
for f in range(1, 5):
    L, x, r0, r1 = letters(f)
    Q = kron_pow(SQX, f)
    D1, D2 = exactD(f, 1), exactD(f, 2)
    P(f"f={f}: |D1 - fib(Q,row0)| = {np.abs(D1 - fib(f, Q, 0)).max():.1e}  |D2 - fib(Q,row1)| = {np.abs(D2 - fib(f, Q, 1)).max():.1e}  |K - D1 D2| check D1D2=D2D1 exact: {np.abs(D1@D2-D2@D1).max():.1e}")
    Y, Yp = rnd(2 ** f), rnd(2 ** f)
    A1, A2 = fib(f, Y, 0), fib(f, Yp, 1)
    P(f"      [A1, A2] = {np.abs(A1 @ A2 - A2 @ A1).max():.1e} (fibrewise, independent Y, Y')")

P("== T1: R8")
for f in range(1, 6):
    L, x, r0, r1 = letters(f)
    n = 4 ** f
    X1 = L["X1"]; pp, pm = (np.eye(n) + X1) / 2, (np.eye(n) - X1) / 2
    models = {"exactQ": (exactD(f, 1), exactD(f, 2)),
              "fibrewise": (fib(f, rnd(2 ** f), 0), fib(f, rnd(2 ** f), 1)),
              "full": (rnd(n), rnd(n))}
    res = []
    for nm, (D1, D2) in models.items():
        # R8 as an operator on x (relative to the size of one half)
        u1 = D1 @ pp
        u2 = D2 @ (pp + pm @ u1)
        u3 = D1 @ (pp @ u2)
        R = pm @ u2 - pm @ u3
        res.append(f"{nm}: |R8|/|pi- u2| = {np.linalg.norm(R) / max(np.linalg.norm(pm @ u2), 1e-300):.2e} (|pi- u2| = {np.linalg.norm(pm@u2)/np.sqrt(n):.2e})")
    P(f"f={f}: " + " | ".join(res))

P("== T2: Proposition S")
for f in (2, 3, 4):
    L, x, r0, r1 = letters(f)
    n = 4 ** f
    X1 = L["X1"]; I = np.eye(n)
    pp, pm = (I + X1) / 2, (I - X1) / 2
    # random perturbations (dense: the identity used only needs pi+- to commute with D1, D2 at Q)
    dc1, dc2, dc3, dB21, dB32, dB31 = (rnd(n) for _ in range(6))
    r0u = [rnd(n) for _ in range(3)]          # arbitrary eps^0 readouts on u1, u2, u3

    def family(D1, D2, s):
        c1 = pp + s * dc1
        u1 = D1 @ c1
        u2 = D2 @ (pp + s * dc2 + (pm + s * dB21) @ u1)
        u3 = D1 @ (s * dc3 + (pp + s * dB32) @ u2 + s * dB31 @ u1)
        return (pm @ u2 - pm @ u3) / s + r0u[0] @ u1 + r0u[1] @ u2 + r0u[2] @ u3

    D1, D2 = exactD(f, 1), exactD(f, 2)
    u1p = D1 @ (pp + pm @ dc1)
    u2p = D2 @ (pp + pm @ dc2 + pm @ u1p + pm @ dB21 @ pp @ u1p)
    u3p = D1 @ (pp @ u2p + pm @ dB32 @ pp @ u2p + pm @ dc3 + pm @ dB31 @ pp @ u1p)
    twin = (pm @ u2p - pm @ u3p) + r0u[0] @ pp @ u1p + r0u[1] @ pp @ u2p + r0u[2] @ pp @ u3p
    line = []
    for s in (1e-2, 1e-4, 1e-6):
        line.append(f"eps={s:g}: |fam-twin|/|twin| = {rel(family(D1, D2, s), twin):.2e}")
    P(f"f={f} exactQ: " + " ; ".join(line))
    F1, F2 = fib(f, rnd(2 ** f), 0), fib(f, rnd(2 ** f), 1)
    P(f"f={f} fibrewise: " + " ; ".join(f"eps={s:g}: |fam|/sqrt(n) = {np.linalg.norm(family(F1, F2, s))/np.sqrt(n):.3e}" for s in (1e-2, 1e-4, 1e-6)))
    G1, G2 = rnd(n), rnd(n)
    P(f"f={f} full: " + " ; ".join(f"eps={s:g}: |fam|/sqrt(n) = {np.linalg.norm(family(G1, G2, s))/np.sqrt(n):.3e}" for s in (1e-2, 1e-4, 1e-6)))
    # extra (not in the claim): the exact twin is a legal program shape 1,2,1; its value is not K I
    K = D1 @ D2
    P(f"f={f}: twin vs K (relative) = {rel(twin, K):.3f} (random perturbations; just shows the twin is a generic exact 3-move output)")

P("== T3: commutant of D1 (fibrewise value)")
for f in range(2, 6):
    L, x, r0, r1 = letters(f)
    n = 4 ** f
    A1 = fib(f, rnd(2 ** f), 0)
    cm = {k: np.abs(A1 @ L[k] - L[k] @ A1).max() for k in ("a1", "b1", "fX2", "X2", "S2", "X1", "a0", "C01")}
    a1, b1 = L["a1"], L["b1"]
    # orders
    def order(g, cap=100000):
        # g is a permutation matrix: order = lcm of cycle lengths
        from math import gcd
        p = np.argmax(np.abs(g), axis=0); seen = np.zeros(n, bool); o = 1
        for s0 in range(n):
            if seen[s0]: continue
            c = 0; t = s0
            while not seen[t]:
                seen[t] = True; t = p[t]; c += 1
            o = o * c // gcd(o, c)
        return o
    P(f"f={f}: [D1fib, .]: " + " ".join(f"{k}={v:.0e}" for k, v in cm.items())
      + f" | ord(a1)={order(a1)} ord(b1)={order(b1)} ord(a1 b1)={order(a1 @ b1)} |[a1,b1]|={np.abs(a1@b1-b1@a1).max():.0f}")
    # lamps
    lamps = [np.linalg.matrix_power(a1, k) @ L["fX2"] @ np.linalg.matrix_power(a1, -k).round() for k in range(min(f, 4))]
    cmax = max(np.abs(p @ q - q @ p).max() for p in lamps for q in lamps)
    P(f"      lamps k<{len(lamps)}: max |[lamp_k, lamp_l]| = {cmax:.0e}, max |lamp^2 - 1| = {max(np.abs(p@p-np.eye(n)).max() for p in lamps):.0e}")

P("== T3b: <L, L'> in GL_f(F2) (L = I + S from a1, L' = I + S^T from b1)")
for f in range(2, 7):
    Lm = np.eye(f, dtype=np.uint8); Lm[np.arange(1, f), np.arange(f - 1)] = 1   # (Lv)_j = v_j + v_(j-1)
    Lp = Lm.T.copy()
    key = lambda M: M.tobytes()
    seen = {key(np.eye(f, dtype=np.uint8))}; frontier = [np.eye(f, dtype=np.uint8)]
    gens = [Lm, Lp]
    while frontier and len(seen) < 3_000_000:
        nf = []
        for M in frontier:
            for g in gens:
                N = (g @ M) % 2
                k = key(N.astype(np.uint8))
                if k not in seen:
                    seen.add(k); nf.append(N.astype(np.uint8))
        frontier = nf
    glsize = 1
    for i in range(f):
        glsize *= (2 ** f - 2 ** i)
    P(f"f={f}: |<L, L'>| = {len(seen)}  |GL_f(F2)| = {glsize}")

P("== T4: braid D_z S_z D_z = S_z D_z S_z")
for f in (1, 2, 3):
    L, x, r0, r1 = letters(f)
    for z in (1, 2, 3):
        D = exactD(f, z); S = L["S%d" % z]
        e_exact = rel(D @ S @ D, S @ D @ S)
        if z == 1:
            Dg = fib(f, rnd(2 ** f), 0)
        elif z == 2:
            Dg = fib(f, rnd(2 ** f), 1)
        else:
            C = L["C01"]; Dg = C @ fib(f, rnd(2 ** f), 0) @ C
        e_gen = rel(Dg @ S @ Dg, S @ Dg @ S)
        P(f"f={f} z={z}: exact {e_exact:.1e}  generic fibrewise {e_gen:.2f}")

P("== T5: Fact A, tau(a1^n) vs 2^(min(2^v,f)-f)")
bad = 0
for f in range(1, 11):
    Lm = np.eye(f, dtype=np.int64); Lm[np.arange(1, f), np.arange(f - 1)] = 1
    for nn in range(1, 25):
        Ln = np.eye(f, dtype=np.int64)
        for _ in range(nn): Ln = (Lm @ Ln) % 2
        # fixed points of v -> Ln v over F2^f
        Mx = (Ln - np.eye(f, dtype=np.int64)) % 2
        # rank over F2
        A = Mx.copy(); r = 0
        for c in range(f):
            piv = [i for i in range(r, f) if A[i, c]]
            if not piv: continue
            A[[r, piv[0]]] = A[[piv[0], r]]
            for i in range(f):
                if i != r and A[i, c]: A[i] ^= A[r]
            r += 1
        fix = f - r
        v = (nn & -nn).bit_length() - 1
        pred = min(2 ** v, f)
        if fix != pred: bad += 1
P(f"Fact A: mismatches over f=1..10, n=1..24: {bad}")

open(__file__.replace(".py", "_out.txt") if len(sys.argv) < 2 else __file__.replace(".py", f"_out_seed{sys.argv[1]}.txt"), "w").write("\n".join(out) + "\n")
