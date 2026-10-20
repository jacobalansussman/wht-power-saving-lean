# b56-compute (wave 16).  CERTIFIED (mod-p rank) demand f(g) = dim(N + gN) - dim N of a coordinate certificate
# N = span{e_t : t in S} in W = Lambda^{a_1} V (x) ... (x) Lambda^{a_r} V, V = C^(2^m), g = K (target) or R = dir(e1).
# Method (b4-cert's identity, re-derived): g^2 is a free shift that keeps N, so N cap gN = (+)_c (N cap W_c), W_c = the
# eigenspaces of g, spanned by Walsh tensors of class c (K: sum of label weights mod 4; R: sum of label bit 0 mod 4).
# Blocks: X-shift character alpha (= xor of all Walsh labels); for even f the Z-shifts Z0(alpha) keep every class inside
# E_alpha, giving characters chi; for odd f, Z_s maps E_alpha -> E_{alpha+s} so f(g) = n * (alpha = 0 part); bit
# permutations (commute with K; those fixing bit 0 commute with R) identify alphas of equal (weight) type.
# Per block: f_block = sum_c rank C[rows of class c] - (#classes - 1) rank C  (rows = Z0-orbit reps of Walsh tuples with
# xor alpha; cols = shift-orbit reps of S with Z0-character chi; C[u,t] = prod det H[U_i, T_i], exact integers).
# CERTIFICATION: every rank is a mod-p rank of (random +-1 compression of) C, so <= the rank over Q; rank C summed over
# blocks is checked to equal dim N exactly, so each rank C is exact and f(g) >= the printed value is PROVED.
# Usage: python3 -I b56.py M A1,A2,... K|R PRED [--comp] [--threads T]
#   PRED: python expression in I (dict region-tuple -> count) e.g. hook "I.get((0,1),0)==1" for (1,k)
#   --comp: use the complement coordinate subspace (same f(g) by duality; pick the smaller)
import sys, os, itertools, time, ctypes
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np
from math import comb
np.seterr(all='ignore')   # Accelerate's matmul raises spurious FP flags; results are checked finite below

HERE = os.path.dirname(os.path.abspath(__file__))

def is_prime(q):
    if q < 2: return False
    d = 2
    while d * d <= q:
        if q % d == 0: return False
        d += 1
    return True

P = (1 << 21) - 1
while not is_prime(P): P -= 2

_lib = None
def lib():
    global _lib
    if _lib is None:
        so = os.path.join(HERE, 'libent.so')
        _lib = ctypes.CDLL(so)
        _lib.fill.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                              ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int]
    return _lib

NTHREADS = 4

def fill(avec, rows, cols):
    """rows: (R, f) uint8 labels, cols: (C, f) uint8 addresses -> (R, C) float64 exact integers"""
    rows = np.ascontiguousarray(rows, dtype=np.uint8); cols = np.ascontiguousarray(cols, dtype=np.uint8)
    out = np.empty((rows.shape[0], cols.shape[0]), dtype=np.float64)
    sz = (ctypes.c_int * len(avec))(*avec)
    if rows.shape[0] and cols.shape[0]:
        lib().fill(len(avec), sz, rows.shape[0], rows.ctypes.data, cols.shape[0], cols.ctypes.data, out.ctypes.data,
                   NTHREADS)
    return out

# ---------------- mod-p rank (blocked, exact in float64: all integers < 2^53) ----------------
def _inv_modp(A, p):
    b = A.shape[0]
    M = np.concatenate([A % p, np.eye(b)], axis=1)
    for c in range(b):
        nz = np.nonzero(M[c:, c])[0]
        assert len(nz), "singular pivot block"
        r = c + nz[0]
        if r != c: M[[c, r]] = M[[r, c]]
        M[c] = (M[c] * pow(int(M[c, c]), p - 2, p)) % p
        col = M[:, c].copy(); col[c] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): M[nzr] = (M[nzr] - np.outer(col[nzr], M[c])) % p
    return M[:, b:]

def _inv_rec(A, p):
    """inverse mod p by recursive 2x2 block (Schur) inversion; valid when all leading principal minors are nonzero
    (true for a pivot block taken in elimination order); falls back to Gauss-Jordan with row swaps otherwise."""
    nA = A.shape[0]
    if nA <= 32:
        return _inv_modp(A, p)
    h = nA // 2
    A1, B, C, D = A[:h, :h], A[:h, h:], A[h:, :h], A[h:, h:]
    try:
        A1i = _inv_rec(A1, p)
    except AssertionError:
        return _inv_modp(A, p)
    CA = (C @ A1i) % p                 # (n-h) x h
    AB = (A1i @ B) % p                 # h x (n-h)
    S = (D - (CA @ B) % p) % p
    try:
        Si = _inv_rec(S, p)
    except AssertionError:
        return _inv_modp(A, p)
    TR = (-(AB @ Si)) % p              # top right
    BL = (-(Si @ CA)) % p              # bottom left
    TL = (A1i - (TR @ CA) % p) % p     # A1i + A1i B Si C A1i = A1i - TR C A1i
    out = np.empty_like(A)
    out[:h, :h] = TL; out[:h, h:] = TR; out[h:, :h] = BL; out[h:, h:] = Si
    return out

def _leaf_pivots(Pn, p):
    """simple elimination on a narrow panel (copy); returns local pivot rows, cols"""
    Pn = Pn.copy()
    piv_rows, piv_cols = [], []
    avail = np.ones(Pn.shape[0], bool)
    for j in range(Pn.shape[1]):
        colj = Pn[:, j]
        cand = np.flatnonzero((colj != 0) & avail)
        if len(cand) == 0: continue
        r = cand[0]; avail[r] = False
        piv_rows.append(r); piv_cols.append(j)
        inv = pow(int(Pn[r, j]), p - 2, p)
        rowv = (Pn[r, j:] * inv) % p
        fac = colj.copy(); fac[r] = 0
        Pn[:, j:] = (Pn[:, j:] - np.outer(fac, rowv)) % p
        Pn[r, j:] = rowv
    return piv_rows, piv_cols

def _elim(M, p, bs):
    """maximal pivot set of M (entries in [0,p)): returns (pivot rows, pivot cols) as indices into M.
    Blocked right-looking elimination; the panel's pivots are found recursively (bs -> 32 -> leaf)."""
    rows = np.arange(M.shape[0]); cols = np.arange(M.shape[1])
    pr_all, pc_all = [], []
    while M.shape[0] and M.shape[1]:
        b = min(bs, M.shape[1])
        if bs <= 32:
            lr, lc = _leaf_pivots(M[:, :b], p)
        else:
            lr, lc = _elim(M[:, :b].copy(), p, 32)
        if len(lr) == 0:
            M = M[:, b:]; cols = cols[b:]; continue
        lr = np.array(lr); lc = np.array(lc)
        pr_all += list(rows[lr]); pc_all += list(cols[lc])
        keep = np.ones(M.shape[0], bool); keep[lr] = False
        oth = np.flatnonzero(keep)
        if M.shape[1] == b or len(oth) == 0: break
        A11 = M[np.ix_(lr, lc)]
        A11i = _inv_rec(A11, p)
        X = (A11i @ M[lr, b:]) % p
        Lf = M[np.ix_(oth, lc)]
        Mn = M[oth, b:]
        del M
        CHR = max(1, int(2e7 // max(1, Mn.shape[1])))
        for st in range(0, Mn.shape[0], CHR):      # chunked Schur update: no full-size temporaries
            T = Lf[st:st + CHR] @ X
            T %= p
            blk = Mn[st:st + CHR]
            blk -= T
            blk %= p
        del Lf
        assert np.isfinite(X).all()
        M = Mn; rows = rows[oth]; cols = cols[b:]
    return pr_all, pc_all

def rank_modp(A, p=P, bsz=256):
    """rank over GF(p), p < 2^21: every float64 sum stays below 256 p^2 < 2^53, so all arithmetic is exact."""
    if isinstance(A, np.ndarray) and A.dtype == np.float64:
        M = A; M %= p            # takes ownership (callers do not reuse A)
    else:
        M = np.array(A, dtype=np.float64) % p
    if M.shape[0] > M.shape[1]: M = np.ascontiguousarray(M.T)
    if M.size == 0: return 0
    pr, pc = _elim(M, p, bsz)
    return len(pr)

def compress_rank(chunks_fn, ncols, nrows, p=P, seed=1):
    """rank lower bound of a tall matrix given by row chunks (generator of float arrays): R C with R random +-1."""
    if nrows <= ncols + 8:
        C = np.concatenate(list(chunks_fn()), axis=0) if nrows else np.zeros((0, ncols))
        return rank_modp(C, p), None
    rng = np.random.default_rng(seed)
    k = ncols + 8
    Y = np.zeros((k, ncols))
    for Cc in chunks_fn():
        Rr = rng.integers(0, 2, size=(k, Cc.shape[0])).astype(np.float64) * 2 - 1
        Y = (Y + Rr @ Cc) % p
        assert np.isfinite(Y).all()
    return rank_modp(Y, p), None

# ---------------- combinatorics on bitmasks ----------------
def popc(x): return bin(x).count('1')

def subsets_masks(n, a):
    out = np.zeros(comb(n, a), dtype=np.uint64)
    for i, T in enumerate(itertools.combinations(range(n), a)):
        v = 0
        for x in T: v |= 1 << x
        out[i] = v
    return out

def xor_perm(masks, z, m):
    """apply x -> x ^ z to the elements of each set (bitmask array)"""
    out = masks.copy()
    for j in range(m):
        if (z >> j) & 1:
            sh = np.uint64(1 << j)
            lo = np.uint64(sum(1 << x for x in range(1 << m) if not (x >> j) & 1))
            out = ((out & lo) << sh) | ((out >> sh) & lo)
    return out

def mask_xor(masks, m):
    """xor of all elements of each set"""
    out = np.zeros(len(masks), dtype=np.int64)
    for x in range(1 << m):
        bit = ((masks >> np.uint64(x)) & np.uint64(1)).astype(bool)
        out[bit] ^= x
    return out

def mask_wsum(masks, m, lab):
    out = np.zeros(len(masks), dtype=np.int64)
    for x in range(1 << m):
        if lab[x] == 0: continue
        bit = ((masks >> np.uint64(x)) & np.uint64(1)).astype(np.int64)
        out += bit * lab[x]
    return out

def masks_to_lists(masks, a, n):
    out = np.zeros((len(masks), a), dtype=np.uint8)
    cnt = np.zeros(len(masks), dtype=np.int64)
    for x in range(n):
        bit = ((masks >> np.uint64(x)) & np.uint64(1)).astype(bool)
        idx = np.nonzero(bit)[0]
        out[idx, cnt[idx]] = x
        cnt[idx] += 1
    assert (cnt == a).all()
    return out

def lexmin_update(best, cand):
    """best, cand: lists of uint64 arrays (one per factor); lexicographic min, in place"""
    better = np.zeros(len(best[0]), bool); eq = np.ones(len(best[0]), bool)
    for b_, c_ in zip(best, cand):
        better |= eq & (c_ < b_); eq &= (c_ == b_)
    for b_, c_ in zip(best, cand): b_[better] = c_[better]

def unique_rows(masklist):
    arr = np.stack(masklist, axis=1)
    _, idx = np.unique(arr, axis=0, return_index=True)
    return np.sort(idx)

# ---------------- the computation ----------------
class Cert:
    def __init__(self, m, avec, pred, comp=False, sgn=None):
        self.m, self.avec, self.r = m, list(avec), len(avec)
        self.sgn = list(sgn) if sgn else [1] * len(avec)
        self.n = 1 << m; self.f = sum(avec); self.pred = pred; self.comp = comp
        n = self.n
        self.subs = {a: subsets_masks(n, a) for a in set(avec)}
        # all tuples as lists of mask arrays (product)
        sizes = [len(self.subs[a]) for a in avec]
        self.tot = int(np.prod(sizes))
        grids = np.meshgrid(*[np.arange(s) for s in sizes], indexing='ij')
        self.all_idx = [g.ravel() for g in grids]
        self.all_masks = [self.subs[a][ix] for a, ix in zip(avec, self.all_idx)]

    def venn_ok(self, masks):
        """evaluate the predicate on every tuple: Venn regions = membership patterns"""
        r = self.r; n = self.n
        regions = [reg for k in range(1, r + 1) for reg in itertools.combinations(range(r), k)]
        cnt = {reg: np.zeros(len(masks[0]), dtype=np.int64) for reg in regions}
        full = np.uint64((1 << n) - 1) if n < 64 else np.uint64(0xFFFFFFFFFFFFFFFF)
        for reg in regions:
            v = np.full(len(masks[0]), full, dtype=np.uint64)
            for i in range(r):
                v = v & (masks[i] if i in reg else (~masks[i] & full))
            # popcount
            c = np.zeros(len(v), dtype=np.int64); w = v.copy()
            while True:
                nzm = w != 0
                if not nzm.any(): break
                c += nzm; w = w & (w - np.uint64(1))
            cnt[reg] = c
        class Lazy(dict):
            pass
        I = {reg: cnt[reg] for reg in regions}
        class Getter:
            def get(self_, k, d=0): return I.get(k, np.zeros(len(masks[0]), dtype=np.int64) + d)
        res = eval(self.pred, {'np': np}, {'I': Getter()})
        res = np.asarray(res, dtype=bool)
        if res.shape == (): res = np.full(len(masks[0]), bool(res))
        return ~res if self.comp else res

    def column_reps(self):
        m, n = self.m, self.n
        ok = self.venn_ok(self.all_masks)
        Sm = [x[ok] for x in self.all_masks]
        self.dimN = int(ok.sum())
        # canonical rep under shifts: lexicographic min over all z (vectorised)
        best = [x.copy() for x in Sm]
        for z in range(1, n):
            lexmin_update(best, [xor_perm(x, z, m) for x in Sm])
        keep = unique_rows(best)
        reps = [b[keep] for b in best]
        self.colmasks = reps
        sig = np.zeros(len(reps[0]), dtype=np.int64)
        for x in reps: sig ^= mask_xor(x, m)
        self.colsig = sig
        self.collist = np.concatenate([masks_to_lists(x, a, n) for x, a in zip(reps, self.avec)], axis=1)
        return len(reps[0])

    def run(self, which='K', verbose=True, maxcols=None):
        m, n, f = self.m, self.n, self.f
        t0 = time.time()
        ncr = self.column_reps()
        if verbose: print(f"  dim N = {self.dimN}, shift-orbit columns {ncr}  ({time.time()-t0:.1f}s)", flush=True)
        # paid element: 'K' = dir(e1)...dir(em); 'R' = dir(e1); 'bIJ..' = dir(e_I)dir(e_J)... (a block of rank len)
        gmask = (n - 1) if which == 'K' else (1 if which == 'R' else sum(1 << int(c) for c in which[1:]))
        lab = [popc(s & gmask) for s in range(n)]
        fs = sum(sg * a for sg, a in zip(self.sgn, self.avec))   # dual factors (Lambda^a V*): g acts as g^-T
        usig = np.zeros(self.tot, dtype=np.int64); ucls = np.zeros(self.tot, dtype=np.int64)
        for x in self.all_masks:
            pass
        for x, sg in zip(self.all_masks, self.sgn):
            usig ^= mask_xor(x, m); ucls += sg * mask_wsum(x, m, lab)
        ucls %= 4
        if which == 'K':
            alphas = [(comb(m, w), (1 << w) - 1) for w in range(m + 1)]
        elif which == 'R':
            alphas = [(comb(m - 1, w), b0 | (((1 << w) - 1) << 1)) for b0 in (0, 1) for w in range(m)]
        else:
            # bit permutations keeping gmask commute with g: classes by (weight inside gmask, weight outside)
            ins = [j for j in range(m) if (gmask >> j) & 1]; out = [j for j in range(m) if not (gmask >> j) & 1]
            alphas = [(comb(len(ins), a) * comb(len(out), b), sum(1 << j for j in ins[:a]) | sum(1 << j for j in out[:b]))
                      for a in range(len(ins) + 1) for b in range(len(out) + 1)]
        if f % 2 == 1:
            alphas = [(n, 0)]
        total = 0; dimchk = 0; maxblk = 0
        for (mult, alpha) in alphas:
            sel = np.nonzero(usig == alpha)[0]
            Z0 = [s for s in range(n) if (fs * popc(s & gmask) + 2 * popc(alpha & s & gmask)) % 4 == 0] if f % 2 == 0 else [0]
            rm = [x[sel] for x in self.all_masks]
            best = [x.copy() for x in rm]
            for s in Z0:
                if s: lexmin_update(best, [xor_perm(x, s, m) for x in rm])
            keep = unique_rows(best)
            rows_m = [x[keep] for x in best]
            rcls = ucls[sel][keep]       # classes are Z0-invariant
            rowlist = np.concatenate([masks_to_lists(x, a, n) for x, a in zip(rows_m, self.avec)], axis=1)
            # Z0 basis for character keys
            basis = []
            span = {0}
            for s in Z0:
                if s not in span:
                    basis.append(s); span = span | {x ^ s for x in span}
            keys = np.zeros(len(self.colsig), dtype=np.int64)
            for i, s in enumerate(basis):
                keys |= ((np.array([popc(int(v) & s) & 1 for v in self.colsig])) << i)
            part = 0
            for kv in np.unique(keys):
                cols = np.nonzero(keys == kv)[0]
                cl = self.collist[cols]
                classes = sorted(set(rcls.tolist()))
                rks = {}
                Ys = []
                colnz = np.zeros(len(cols), bool)
                for c in classes:
                    ridx = np.nonzero(rcls == c)[0]
                    rl = rowlist[ridx]
                    CH = max(1, int(4e6 // max(1, len(cols))))
                    def gen(rl=rl, CH=CH):
                        for st in range(0, len(rl), CH):
                            Cc = fill(self.avec, rl[st:st + CH], cl)
                            colnz[:] |= (Cc != 0).any(axis=0)
                            yield Cc
                    rk, Y = compress_rank(gen, len(cols), len(rl), seed=1 + c)
                    rks[c] = rk; del Y
                    maxblk = max(maxblk, len(cols))
                # rank C = number of nonzero columns EXACTLY: the columns are projections of distinct shift orbits,
                # with disjoint supports in the standard basis [proved]; so no rank computation is needed.
                rkC = int(colnz.sum())
                if len(classes) == 2:
                    blk = rks[classes[0]] + rks[classes[1]] - rkC
                elif len(classes) == 1:
                    blk = 0          # the whole block lies in one eigenspace of g: N_block is g-stable
                else:
                    raise RuntimeError("more than two classes in one X-character block: impossible (g^2 is a scalar there)")
                part += blk; dimchk += mult * rkC
            total += mult * part
            if verbose:
                print(f"  alpha={alpha:0{m}b} x{mult}: |Z0|={len(Z0)}, rows {len(keep)}, part {part}  ({time.time()-t0:.1f}s)",
                      flush=True)
        return self.dimN, dimchk, total, maxblk

if __name__ == '__main__':
    args = sys.argv[1:]
    comp = '--comp' in args
    if '--threads' in args:
        NTHREADS = int(args[args.index('--threads') + 1])
    args = [a for i, a in enumerate(args) if not a.startswith('--') and not (i > 0 and args[i - 1] == '--threads')]
    m = int(args[0]); which = args[2]; pred = args[3]
    toks = args[1].split(',')     # 'd3' = Lambda^3 V* (dual factor: letters act by inverse transpose)
    avec = [int(x.lstrip('d')) for x in toks]; sgn = [-1 if x.startswith('d') else 1 for x in toks]
    t0 = time.time()
    print(f"m={m} W=Lambda^{toks} g={which} N: {pred}{' (complement)' if comp else ''}  p={P}", flush=True)
    C = Cert(m, avec, pred, comp, sgn)
    dimN, dimchk, val, mb = C.run(which)
    ok = (dimchk == dimN)
    print(f"RESULT m={m} {avec} {which} [{pred}]{' comp' if comp else ''}: dim N {dimN}, rank check {dimchk} "
          f"{'EXACT' if ok else 'MISMATCH'}, f({which}) >= {val}{' (certified)' if ok else ' (NOT certified)'}"
          f"  max block cols {mb}  ({time.time()-t0:.1f}s)", flush=True)
