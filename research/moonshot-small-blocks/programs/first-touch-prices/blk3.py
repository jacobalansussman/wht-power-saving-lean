# blk3.py -- w7-mixed3 (15 Oct 2026).  Block-3 mixed programs: one spec, three readings.
#   spec = (W, targets, moves); targets[w] in {'K', 'R'}; a move = (zdir, items), items = [(coef, letter, src)],
#   src < W is the raw input x_src, src >= W is generator src (moves are numbered W, W+1, ...).
#   The move's generator is u = dir(zdir) * sum coef * letter * gen[src].
# Readings:
#   rep_gens(spec, f, coefs)  -> per generator, a list of W matrices rho_f(slot coefficient)  (complex floats)
#   f2test(...)               -> the f = 2 widened row-space test (Sym^2 = D8 + S28, Lambda^2 = A28; each irreducible with
#                                multiplicity 1 under Z4 wr S8, so the left C[B']-span is a row space per block)
#   formal(spec, p)           -> mixlib.Program (exact group elements, coefficients mod p) for hook potentials / exact solve
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

# ---------------------------------------------------------------- exact-ish 8x8 matrices as complex floats ----------
def dirm(z):
    M = np.zeros((8, 8), complex)
    for x in range(8):
        M[x, x] += (1 + 1j) / 2; M[x ^ z, x] += (1 - 1j) / 2
    return M
D = {1: dirm(1), 2: dirm(2), 4: dirm(4)}
KM = D[1] @ D[2] @ D[4]
RM = D[1] @ D[2]
TGT = {'K': KM, 'R': RM}

def mono(perm, ph):
    M = np.zeros((8, 8), complex)
    for x in range(8): M[perm[x], x] = 1j ** (ph[x] % 4)
    return M

def is_mono_z4(M, tol=1e-9):
    nz = np.abs(M) > tol
    if not ((nz.sum(0) == 1).all() and (nz.sum(1) == 1).all()): return False
    v = M[nz]
    return bool(np.all(np.min(np.abs(v[:, None] - np.array([1, 1j, -1, -1j])[None, :]), axis=1) < tol))

def rand_widened(rng):
    perm = list(rng.permutation(8)); ph = list(rng.integers(0, 4, 8))
    return mono(perm, ph)

def affine_perm(A, c):
    out = []
    for x in range(8):
        y = c
        for j in range(3):
            if x >> j & 1: y ^= A[j]
        out.append(y)
    return out

def rand_engine(rng):
    while True:
        A = [int(rng.integers(8)) for _ in range(3)]
        if len({A[0], A[1], A[2], A[0] ^ A[1], A[0] ^ A[2], A[1] ^ A[2], A[0] ^ A[1] ^ A[2]} - {0}) == 7: break
    perm = affine_perm(A, int(rng.integers(8)))
    a = rng.integers(0, 4, 3); b = rng.integers(0, 2, 3); c0 = int(rng.integers(4))
    ph = []
    for x in range(8):
        xs = [(x >> j) & 1 for j in range(3)]
        v = c0 + sum(int(a[j]) * xs[j] for j in range(3)) + 2 * (b[0] * xs[0] * xs[1] + b[1] * xs[0] * xs[2] + b[2] * xs[1] * xs[2])
        ph.append(int(v) % 4)
    return mono(perm, ph)

I8 = np.eye(8, dtype=complex)

# ---------------------------------------------------------------- f = 2 blocks ---------------------------------------
def blocks8():
    Dg = np.zeros((8, 64)); S = []; A = []
    for x in range(8): Dg[x, x * 8 + x] = 1
    for x in range(8):
        for y in range(x + 1, 8):
            v = np.zeros(64); v[x * 8 + y] = v[y * 8 + x] = 1 / np.sqrt(2); S.append(v)
            w = np.zeros(64); w[x * 8 + y] = 1 / np.sqrt(2); w[y * 8 + x] = -1 / np.sqrt(2); A.append(w)
    return {'D': Dg.astype(complex), 'S': np.array(S, complex), 'A': np.array(A, complex)}
QB = blocks8()

def kronf(M, f):
    out = M
    for _ in range(f - 1): out = np.kron(out, M)
    return out

# ---------------------------------------------------------------- readings -------------------------------------------
def rep_gens(spec, f, coefs=None):
    """coefs: optional flat complex vector overriding every item's coef (in order)"""
    W, targets, moves = spec
    n = 8 ** f
    gens = [[(np.eye(n, dtype=complex) if l == w else None) for l in range(W)] for w in range(W)]
    Df = {z: kronf(D[z], f) for z in (1, 2, 4)}
    k = 0
    for (z, items) in moves:
        acc = [None] * W
        for (c, L, src) in items:
            if coefs is not None: c = coefs[k]
            k += 1
            Lf = kronf(L, f)
            for l in range(W):
                g = gens[src][l]
                if g is None: continue
                t = c * (Lf @ g)
                acc[l] = t if acc[l] is None else acc[l] + t
        gens.append([None if a is None else Df[z] @ a for a in acc])
    return gens

def ncoefs(spec):
    return sum(len(items) for (_, items) in spec[2])

def f2_residual(spec, coefs=None, tol=1e-9, detail=False, gens=None):
    W, targets, moves = spec
    if gens is None: gens = rep_gens(spec, 2, coefs)
    tot = 0.0; out = {}
    for w in range(W):
        T2 = np.kron(TGT[targets[w]], TGT[targets[w]])
        for b, Q in QB.items():
            rows = []
            for g in gens:
                blk = np.hstack([Q @ (g[l] if g[l] is not None else np.zeros((64, 64))) for l in range(W)])
                rows.append(blk)
            S = np.vstack(rows)
            T = np.zeros((Q.shape[0], 64 * W), complex); T[:, 64 * w:64 * (w + 1)] = Q @ T2
            U, s, Vh = np.linalg.svd(S, full_matrices=False)
            r = int((s > tol * s[0]).sum()) if s.size else 0
            Vr = Vh[:r]
            Rr = T - (T @ Vr.conj().T) @ Vr
            res = float(np.sum(np.abs(Rr) ** 2) / np.sum(np.abs(T) ** 2))
            out[(w, b)] = (res, r); tot += res
    return (tot, out) if detail else tot

def f2_lsq_residual_vec(spec, coefs, wslots=None):
    """smooth residual vector: for each target slot and block, T - X S with X the least-squares readout (used by LM)"""
    W, targets, moves = spec
    gens = rep_gens(spec, 2, coefs)
    vecs = []
    for w in (range(W) if wslots is None else wslots):
        T2 = np.kron(TGT[targets[w]], TGT[targets[w]])
        for b, Q in QB.items():
            S = np.vstack([np.hstack([Q @ (g[l] if g[l] is not None else np.zeros((64, 64))) for l in range(W)]) for g in gens])
            T = np.zeros((Q.shape[0], 64 * W), complex); T[:, 64 * w:64 * (w + 1)] = Q @ T2
            # least squares X S = T  -> X = T S^+ (ridge for stability)
            St = S.conj().T
            G = S @ St + 1e-10 * np.eye(S.shape[0])
            X = np.linalg.solve(G.T, (T @ St).T).T
            Rr = (T - X @ S) / np.sqrt(np.sum(np.abs(T) ** 2))
            vecs.append(Rr.ravel())
    v = np.concatenate(vecs)
    return np.concatenate([v.real, v.imag])

# ---------------------------------------------------------------- formal (exact) reading ----------------------------
def to_gm(M):
    """complex 8x8 matrix with entries in Z[i]/2^s -> hooklib.GM (exact)"""
    from hooklib import GM
    for s in range(0, 8):
        A = M * (2 ** s)
        re = np.rint(A.real); im = np.rint(A.imag)
        if np.allclose(re, A.real, atol=1e-9) and np.allclose(im, A.imag, atol=1e-9):
            return GM(re.astype(np.int64), im.astype(np.int64), s)
    raise ValueError('not dyadic')

def formal(spec, p, coef_modp=None):
    """exact program over F_p: coefficients must be given mod p (coef_modp: list of ints) or be small Gaussian integers"""
    import mixlib
    W, targets, moves = spec
    prog = mixlib.Program(W, p)
    F = mixlib.Field(p)
    k = 0
    for (z, items) in moves:
        its = []
        for (c, L, src) in items:
            if coef_modp is not None: cm = coef_modp[k]
            else:
                cr, ci = int(round(c.real)), int(round(c.imag))
                assert abs(c - (cr + 1j * ci)) < 1e-9, c
                cm = (cr + ci * F.I) % p
            k += 1
            its.append((cm, to_gm(L), src))
        prog.move(to_gm(D[z]), its)
    return prog

def targets_gm(spec):
    from hooklib import gm_dir
    W, targets, moves = spec
    D1, D2, D4 = gm_dir(1), gm_dir(2), gm_dir(4)
    T = {'K': D1 @ D2 @ D4, 'R': D1 @ D2}
    return {w: T[targets[w]] for w in range(W)}

# ---------------------------------------------------------------- plain programs (controls) --------------------------
def plain_spec(nK, nR, drop_last=False):
    W = nK + nR
    targets = ['K'] * nK + ['R'] * nR
    moves = []
    nxt = W
    for w in range(nK):
        moves.append((1, [(1, I8, w)])); a = nxt; nxt += 1
        moves.append((2, [(1, I8, a)])); b = nxt; nxt += 1
        moves.append((4, [(1, I8, b)])); nxt += 1
    for w in range(nK, W):
        moves.append((1, [(1, I8, w)])); a = nxt; nxt += 1
        if drop_last and w == W - 1: break
        moves.append((2, [(1, I8, a)])); nxt += 1
    return (W, targets, moves)

if __name__ == '__main__':
    for nK, nR, dl in [(1, 5, False), (1, 5, True), (2, 2, False), (2, 2, True)]:
        sp = plain_spec(nK, nR, dl)
        tot, out = f2_residual(sp, detail=True)
        print(f'plain {nK}+{nR} drop_last={dl}: moves {len(sp[2])}, f2 residual {tot:.3e}',
              {k: round(v[0], 4) for k, v in out.items() if v[0] > 1e-8})
