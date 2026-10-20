"""h-L1: exact numerical check of the slot-residual lemma (L1) in the semantics of our engine, before formalising it.

L1 (statement checked here).  For a stage-B frame map X : BXMap H alpha of the window kind used by the #401/#588 network
(X = BXMap.ofConj d q: alpha = H (+) beta, liftG G = conjPerm d G : y |-> inlPerm G (d^T y), liftS s = s (+) q) and an
entrance label U = (G, s) of H (dim U = |s| = e):

    X.Phi U.top = C * blockMat z * C' * X.Phi U,     C = kg G',  C' = kg (apInv G'),  G' = X.liftG G,
                                                    z = the unit vectors e_{inl k}, k not in s   (rank |H| - e),
    C' * C = 1 = C * C'.

Engine semantics (all from the Lean definitions, OAI BinaryFrames.lean and Work/GFrame/Clifford/Core.lean):
  Space alpha = F_2^alpha (ints), dot = parity of AND, sign t = (-1)^t, tint t = i^t, lum x = prod tint(x_k) = i^|x|,
  wal x y = sign(dot x y), walI = digitProd littleH' = wal / 2^n, wrap f = wal * diag f * walI,
  pk S = wrap (x |-> prod_{k in S} tint (x k)),  kg G x y = [G x = y] lum x / lum y,  core G S = kg G * pk S * kg G^-1,
  shift z i j = [i + z = j],  dir z = ca*1 + cb*shift z,  ca = (1+i)/2, cb = (1-i)/2,  blockMat z = prod_i dir (z i).
Arithmetic: exact Gaussian integers (int64 real and imaginary parts) with explicit powers of 2 (2*ca = 1+i,
2*cb = 1-i, 2^n * walI = wal).  No floating point in the main checks; small sizes are cross-checked against dense
matrix products of the definitions (complex128; every entry there is a dyadic number with < 30 significant bits, so the
products are exact, and we additionally compare after scaling to integers).

Shapes: |H| = 9 (the #588 / U8 label space) with entrance dimensions e = 0,1,2,4,5 (slot widths 9,8,7,5,4: all widths
that occur in #588's 616 slots and in the U8 word), plus every e in 0..9; small replicas |H| = 2..4 with every e and
dense cross-checks.  Controls that MUST fail are run on every instance (wrong directions, missing adapters, swapped
adapters, a different window's adapter).
"""
import sys, json, random, time
sys.path.append('<verify>/checks/wht33/tinygame/pylib')   # numpy 2.0.2 (local copy)
import numpy as np

# ---------------------------------------------------------------- exact Gaussian-integer matrices ------------------
class GI:
    __slots__ = ('re', 'im')
    def __init__(self, re, im=None):
        self.re = np.asarray(re, dtype=np.int64)
        self.im = np.zeros_like(self.re) if im is None else np.asarray(im, dtype=np.int64)
    def copy(self): return GI(self.re.copy(), self.im.copy())
    def __add__(a, b): return GI(a.re + b.re, a.im + b.im)
    def __sub__(a, b): return GI(a.re - b.re, a.im - b.im)
    def scal(a, k):               # multiply by the integer k
        return GI(a.re * k, a.im * k)
    def mul_1pi(a):               # (1 + i) * a
        return GI(a.re - a.im, a.re + a.im)
    def mul_1mi(a):               # (1 - i) * a
        return GI(a.re + a.im, a.im - a.re)
    def eq(a, b): return bool(np.array_equal(a.re, b.re) and np.array_equal(a.im, b.im))
    def maxabs(a): return int(max(np.abs(a.re).max(), np.abs(a.im).max()))

def unit_pow(k):
    """i^k for an int array k -> (re, im) int arrays"""
    k = np.asarray(k) % 4
    re = np.where(k == 0, 1, np.where(k == 2, -1, 0)).astype(np.int64)
    im = np.where(k == 1, 1, np.where(k == 3, -1, 0)).astype(np.int64)
    return re, im

def rows_times_units(M, k):       # row x multiplied by i^{k[x]}
    ur, ui = unit_pow(k)
    ur = ur[:, None]; ui = ui[:, None]
    return GI(M.re * ur - M.im * ui, M.re * ui + M.im * ur)

def cols_times_units(M, k):       # column y multiplied by i^{k[y]}
    ur, ui = unit_pow(k)
    ur = ur[None, :]; ui = ui[None, :]
    return GI(M.re * ur - M.im * ui, M.re * ui + M.im * ur)

# ---------------------------------------------------------------- the space F_2^n --------------------------------
def popcount_table(N):
    pc = np.zeros(N, dtype=np.int64)
    for x in range(1, N):
        pc[x] = pc[x >> 1] + (x & 1)
    return pc

def linear_table(images, n):
    """table of the F_2-linear map sending unit vector e_j to images[j] (ints), on all 2^n vectors"""
    N = 1 << n
    T = np.zeros(N, dtype=np.int64)
    for y in range(1, N):
        low = (y & -y).bit_length() - 1
        T[y] = T[y & (y - 1)] ^ images[low]
    return T

def is_bijection(T):
    return len(np.unique(T)) == len(T)

# ---------------------------------------------------------------- engine objects ----------------------------------
class Space:
    def __init__(self, n):
        self.n = n; self.N = 1 << n
        self.pc = popcount_table(self.N)
        self.idx = np.arange(self.N, dtype=np.int64)

    # kg G as a LEFT factor:  (kg G * M)(x, .) = lum x / lum (G x) * M (G x, .)        [kg_mul_apply]
    def kg_left(self, Gt, M):
        R = GI(M.re[Gt, :], M.im[Gt, :])
        return rows_times_units(R, self.pc - self.pc[Gt])

    # kg G as a RIGHT factor: (M * kg G)(., y) = M(., G^-1 y) * lum (G^-1 y) / lum y
    def kg_right(self, M, Gt):
        Ginv = np.empty_like(Gt); Ginv[Gt] = self.idx
        R = GI(M.re[:, Ginv], M.im[:, Ginv])
        return cols_times_units(R, self.pc[Ginv] - self.pc)

    # 2^n * pk S, from pk S = wal * diag f * walI = (g(x xor y))/2^n, g = Walsh transform of f, f(u) = i^{|u & S|}
    def pk_scaled(self, S):
        fr, fi = unit_pow(self.pc[self.idx & S])
        g_re, g_im = fr.copy(), fi.copy()
        h = 1
        while h < self.N:                       # in-place fast Walsh-Hadamard transform (exact, integers)
            for arr in (g_re, g_im):
                a = arr.reshape(-1, 2 * h)
                x = a[:, :h].copy(); y = a[:, h:].copy()
                a[:, :h] = x + y; a[:, h:] = x - y
            h *= 2
        X = self.idx[:, None] ^ self.idx[None, :]
        return GI(g_re[X], g_im[X])

    # 2^n * core G S = kg G * (2^n pk S) * kg (apInv G)
    def core_scaled(self, Gt, S):
        Ginv = np.empty_like(Gt); Ginv[Gt] = self.idx
        return self.kg_right(self.kg_left(Gt, self.pk_scaled(S)), Ginv)

    # (2 dir z) * M = (1+i) M + (1-i) shift z * M,  (shift z * M)(x, .) = M(x + z, .)
    def dir2_left(self, z, M):
        Sh = GI(M.re[self.idx ^ z, :], M.im[self.idx ^ z, :])
        return M.mul_1pi() + Sh.mul_1mi()

    def block2_left(self, zs, M):   # 2^r blockMat zs * M  (List.ofFn prod: dir z0 * dir z1 * ...; applied right to left)
        for z in reversed(zs):
            M = self.dir2_left(z, M)
        return M

    def eye(self):
        return GI(np.eye(self.N, dtype=np.int64))

# ---------------------------------------------------------------- dense definitional cross-check (small n) --------
def dense_defs(sp, Gt, S):
    """core G S by literal matrix products of the definitions (complex128), for small n"""
    N, n = sp.N, sp.n
    I = 1j
    lum = np.array([I ** int(sp.pc[x]) for x in range(N)])
    wal = np.array([[(-1) ** int(sp.pc[x & y]) for y in range(N)] for x in range(N)], dtype=complex)
    walI = wal / 2 ** n
    f = np.array([I ** int(sp.pc[u & S]) for u in range(N)])
    pk = wal @ np.diag(f) @ walI
    kg = np.zeros((N, N), dtype=complex)
    for x in range(N):
        kg[x, Gt[x]] = lum[x] / lum[Gt[x]]
    Ginv = np.empty_like(Gt); Ginv[Gt] = sp.idx
    kgi = np.zeros((N, N), dtype=complex)
    for x in range(N):
        kgi[x, Ginv[x]] = lum[x] / lum[Ginv[x]]
    return kg @ pk @ kgi, kg, kgi

def dense_dir(sp, z):
    ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
    Sh = np.zeros((sp.N, sp.N), dtype=complex)
    for i in range(sp.N):
        Sh[i, i ^ z] = 1
    return ca * np.eye(sp.N) + cb * Sh

def to_complex(M, scale):
    return (M.re + 1j * M.im) / float(scale)

# ---------------------------------------------------------------- random windows and labels ----------------------
def rand_invertible(h, rng):
    while True:
        cols = [rng.randrange(1, 1 << h) for _ in range(h)]
        if is_bijection(linear_table(cols, h)):
            return cols

def rand_orthogonal(n, rng, steps=None):
    """random F_2-orthogonal n x n matrix (d^T d = 1): product of coordinate permutations and transvections
    x |-> x + (v.x) v with v of even weight; returned as a list of rows (ints, bit j of row i = d[i][j])"""
    d = [1 << i for i in range(n)]
    def mul(A, B):            # A * B, rows as ints
        out = []
        for i in range(n):
            r = 0
            for j in range(n):
                if A[i] >> j & 1:
                    r ^= B[j]
            out.append(r)
        return out
    for _ in range(steps or 4 * n):
        if rng.random() < 0.4:
            p = list(range(n)); rng.shuffle(p)
            P = [1 << p[i] for i in range(n)]
            d = mul(P, d)
        else:
            while True:
                v = rng.randrange(1, 1 << n)
                if bin(v).count('1') % 2 == 0:
                    break
            T = [(1 << i) ^ (v if v >> i & 1 else 0) for i in range(n)]
            d = mul(T, d)
    # check orthogonality d^T d = 1  <=>  rows orthonormal for the standard dot
    for i in range(n):
        for j in range(n):
            assert bin(d[i] & d[j]).count('1') % 2 == (1 if i == j else 0)
    return d

def lift_table(h, b, d, Gcols):
    """table of liftG G = conjPerm d G : y |-> inlPerm G (d^T y) on F_2^(h+b); bits 0..h-1 = inl, h.. = inr"""
    n = h + b
    # d^T y : (d^T y)_i = sum_j d[j][i] y_j  -> image of e_j is row j of d
    dT_images = [d[j] for j in range(n)]
    DT = linear_table(dT_images, n)
    Gt_h = linear_table(Gcols, h)
    mask_h = (1 << h) - 1
    return (Gt_h[DT & mask_h]) | (DT & ~mask_h)

# ---------------------------------------------------------------- one instance ---------------------------------------
def check_instance(h, b, e, rng, dense=False):
    n = h + b; sp = Space(n)
    d = rand_orthogonal(n, rng)
    while True:                    # a non-identity window (so the adapters matter)
        if d != [1 << i for i in range(n)]:
            break
        d = rand_orthogonal(n, rng)
    Gcols = rand_invertible(h, rng)
    s = sorted(rng.sample(range(h), e))
    q = [j for j in range(b) if rng.random() < 0.5]
    smask = sum(1 << k for k in s)
    qmask = sum(1 << (h + j) for j in q)
    T_s = smask | qmask                               # liftS s = s (+) q
    T_u = ((1 << h) - 1) | qmask                      # liftS univ
    missing = T_u & ~T_s
    assert missing == ((1 << h) - 1) & ~smask        # = inl '' (univ \ s): the images of the missing label coordinates
    zs = [1 << k for k in range(n) if missing >> k & 1]
    r = len(zs); assert r == h - e
    Gt = lift_table(h, b, d, Gcols); assert is_bijection(Gt)
    Ginv = np.empty_like(Gt); Ginv[Gt] = sp.idx

    core_s = sp.core_scaled(Gt, T_s)                  # 2^n Phi U
    core_u = sp.core_scaled(Gt, T_u)                  # 2^n Phi U.top
    LHS = core_u.scal(1 << r)

    def rhs(C_tab, Ci_tab, dirs):
        M = sp.kg_left(Ci_tab, core_s)                # C' * Phi U
        M = sp.block2_left(dirs, M)                   # 2^r blockMat z * ...
        return sp.kg_left(C_tab, M)                   # C * ...
    res = {}
    res['L1'] = LHS.eq(rhs(Gt, Ginv, zs))
    # C' * C = 1 and C * C' = 1
    E = sp.eye()
    res['CpC_one'] = sp.kg_left(Ginv, sp.kg_left(Gt, E)).eq(E)
    res['CCp_one'] = sp.kg_left(Gt, sp.kg_left(Ginv, E)).eq(E)
    # pk (missing) = blockMat z   (pk_eq_blockMat): 2^n pk = 2^(n-r) * (2^r blockMat z)
    res['pk_eq_block'] = sp.pk_scaled(missing).eq(sp.block2_left(zs, E).scal(1 << (n - r)))
    # controls: each must FAIL (when the control is applicable)
    ctrl = {}
    if e > 0:
        k = rng.choice(s); bad = zs[:-1] + [1 << k]
        ctrl['dir_already_applied'] = not LHS.eq(rhs(Gt, Ginv, bad))
    if b > 0:
        j = rng.randrange(b); bad = zs[:-1] + [1 << (h + j)]
        ctrl['dir_extra_coordinate'] = not LHS.eq(rhs(Gt, Ginv, bad))
    if r >= 1:                                    # with r = 0 the residual is the identity: adapter controls N/A
        Id = sp.idx.copy()
        ctrl['no_adapters'] = not LHS.eq(rhs(Id, Id, zs))
        ctrl['swapped_adapters'] = not LHS.eq(rhs(Ginv, Gt, zs))
        d2 = rand_orthogonal(n, rng)
        Gt2 = lift_table(h, b, d2, Gcols); Gi2 = np.empty_like(Gt2); Gi2[Gt2] = sp.idx
        if not np.array_equal(Gt2, Gt):
            ctrl['other_window_adapter'] = not LHS.eq(rhs(Gt2, Gi2, zs))
    if r >= 1:
        ctrl['rank_one_short'] = not LHS.eq(rhs(Gt, Ginv, zs[:-1]).scal(2))
    res['controls_fail_as_expected'] = all(ctrl.values())
    res['controls'] = ctrl
    if dense:
        with np.errstate(all='ignore'):   # numpy 2.0 + Accelerate raise spurious FP flags in complex matmul; compared exactly
            Pu, kg, kgi = dense_defs(sp, Gt, T_u)
            Ps, _, _ = dense_defs(sp, Gt, T_s)
            B = np.eye(sp.N, dtype=complex)
            for z in zs:
                B = B @ dense_dir(sp, z)
            res['dense_core_u'] = bool(np.array_equal(Pu * 2 ** n, to_complex(core_u, 1)))
            res['dense_core_s'] = bool(np.array_equal(Ps * 2 ** n, to_complex(core_s, 1)))
            res['dense_L1'] = bool(np.array_equal(Pu * 2 ** (n + r), (kg @ B @ kgi @ Ps) * 2 ** (n + r)))
            res['dense_CpC'] = bool(np.array_equal(kgi @ kg, np.eye(sp.N)))
    res['shape'] = dict(h=h, b=b, n=n, e=e, r=r, q=q, s=s)
    res['max_entry'] = LHS.maxabs()
    return res

def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 20261013)
    big_b = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1]
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    out = {'instances': [], 'summary': {}}
    t0 = time.time()
    plan = []
    for h in (2, 3, 4):                          # small replicas, every entrance dimension, dense cross-checks
        for b in (1, 2):
            for e in range(h + 1):
                for _ in range(2):
                    plan.append((h, b, e, True))
    for b in big_b:                               # the real label space |H| = 9
        for e in range(10):
            for _ in range(reps if e in (0, 1, 2, 4, 5) else 1):
                plan.append((9, b, e, False))
    ok = True; fails = []; nctrl = [0, 0]; coinc = []
    for (h, b, e, dense) in plan:
        r = check_instance(h, b, e, rng, dense)
        good = r['L1'] and r['CpC_one'] and r['CCp_one'] and r['pk_eq_block']
        if h == 9:                                # at the real size every control must fail
            good = good and r['controls_fail_as_expected']
        nctrl[0] += len(r['controls']); nctrl[1] += sum(1 for v in r['controls'].values() if v)
        for k, v in r['controls'].items():
            if not v:
                coinc.append(dict(r['shape'], control=k))
        if dense:
            good = good and r['dense_core_u'] and r['dense_core_s'] and r['dense_L1'] and r['dense_CpC']
        ok &= good
        if not good:
            fails.append(r)
        out['instances'].append(r)
        print(f"h={h} b={b} e={e} r={h-e} dense={dense}: L1={r['L1']} C'C=1:{r['CpC_one']} pk=block:{r['pk_eq_block']} "
              f"controls={r['controls']} {'dense ok' if dense and r['dense_L1'] and r['dense_core_u'] else ''}", flush=True)
    occ = [x for x in out['instances'] if x['shape']['h'] == 9 and x['shape']['e'] in (0, 1, 2, 4, 5)]
    out['summary'] = dict(all_pass=ok, n_instances=len(out['instances']), n_fail=len(fails),
                          h9_occurring_shapes=len(occ), h9_occurring_pass=sum(1 for x in occ if x['L1']),
                          controls_run=nctrl[0], controls_failed_as_expected=nctrl[1],
                          controls_coincidentally_true=coinc, seconds=round(time.time() - t0, 1))
    print(json.dumps(out['summary']))
    return out, ok

if __name__ == '__main__':
    out, ok = main()
    path = sys.argv[4] if len(sys.argv) > 4 else None
    if path:
        json.dump(out, open(path, 'w'), indent=1, default=str)
    print('ALL L1 CHECKS PASS' if ok else 'SOME L1 CHECK FAILED')
    sys.exit(0 if ok else 1)
