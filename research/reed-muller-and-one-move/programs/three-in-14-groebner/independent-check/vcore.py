"""vcore.py (v-r3-3in14, skeptic): my own code, written from the definitions in rm-theory 8.5-8.7; no code shared
with r3-3in14, rm2-solve or v-rm2-solve.

Field GF(32) = F2[x]/(x^5+x^2+1), alpha = x (the same field and primitive element as everyone, so that target labels
k <-> c = alpha^(2k) can be compared).  Everything else is basis-free:
- a symmetric bilinear form on GF(32) is stored by its values on the POLYNOMIAL basis e_i = alpha^i (not a self-dual
  basis), as a 15-bit vector over the monomials (i, j), i <= j;
- a direction x (nonzero field element) gives the form b_x(u, v) = Tr(xu) Tr(xv)  (= z z^T in any self-dual basis);
- target c gives the form Tr(c u v)  (= T_k in a self-dual basis when c = alpha^(2k));
- a "form" Q is any of the 2^15 linear functionals on that space; Q(x) := Q(b_x), p_c(Q) := Q(Tr(c..)).
The generic tests below take a list of target vectors and move vectors in ANY 15-bit coordinates of Sym (so they also
run on the other tracks' instances given in their own self-dual coordinates, using their own matrices).
"""
import itertools
import numpy as np

M5 = 0b100101
EXP = [0] * 62; LOG = {}
_x = 1
for _k in range(31):
    EXP[_k] = _x; LOG[_x] = _k
    _x <<= 1
    if _x & 32: _x ^= M5
assert _x == 1
for _k in range(31, 62): EXP[_k] = EXP[_k - 31]


def mul(a, b):
    return 0 if a == 0 or b == 0 else EXP[LOG[a] + LOG[b]]


def tr(a):
    s = 0; y = a
    for _ in range(5):
        s ^= y; y = mul(y, y)
    assert s in (0, 1)
    return s


MONS = [(i, j) for i in range(5) for j in range(i, 5)]      # 15 monomials, my own order
POLYB = [EXP[i] for i in range(5)]


def sym_vec_from_bil(f):
    """f(u, v) -> 15-bit vector (values on the polynomial basis)."""
    v = 0
    for t, (i, j) in enumerate(MONS):
        if f(POLYB[i], POLYB[j]): v |= 1 << t
    return v


def dir_vec(x):
    return sym_vec_from_bil(lambda u, v: tr(mul(x, u)) & tr(mul(x, v)))


def tgt_vec(c):
    return sym_vec_from_bil(lambda u, v: tr(mul(c, mul(u, v))))


def mat_to_vec(T):
    """symmetric 0/1 matrix in some coordinates -> 15-bit vector in the same coordinates (monomials i <= j)."""
    v = 0
    for t, (i, j) in enumerate(MONS):
        if T[i][j] & 1: v |= 1 << t
    return v


def zz_vec(z):
    """move direction given as a 5-bit coordinate vector z (bit j = coordinate j) -> vector of z z^T."""
    bits = [(z >> j) & 1 for j in range(5)]
    return mat_to_vec([[bits[i] & bits[j] for j in range(5)] for i in range(5)])


LAM = np.arange(1 << 15, dtype=np.int64)
_PAR = np.zeros(1 << 15, dtype=np.uint8)
for _b in range(15):
    _PAR ^= ((LAM >> _b) & 1).astype(np.uint8)


def parity_table(vec):
    """array over all 2^15 functionals Q of Q(vec) = parity(Q & vec)."""
    return _PAR[LAM & vec]


class Inst:
    """W targets and M moves, all as 15-bit vectors in one coordinate system."""
    def __init__(self, tvecs, mvecs):
        self.W = len(tvecs); self.M = len(mvecs)
        self.tv = list(tvecs); self.mv = list(mvecs)
        self.Qm = np.stack([parity_table(v) for v in mvecs])            # M x 2^15
        self.Pt = np.stack([parity_table(v) for v in tvecs])            # W x 2^15
        self.cnt = self.Qm.sum(0).astype(np.int32)
        self.r = self.Pt.sum(0).astype(np.int32)
        self.slack = self.cnt - self.r
        # move-support bitmask of each form, parity bitmask of each form
        w = (np.int64(1) << np.arange(self.M, dtype=np.int64))
        self.Dmask = (self.Qm.astype(np.int64) * w[:, None]).sum(0)
        wl = (np.int64(1) << np.arange(self.W, dtype=np.int64))
        self.Pmask = (self.Pt.astype(np.int64) * wl[:, None]).sum(0)

    def deletion_ok(self):
        return bool((self.slack >= 0).all())

    def tight_supports(self):
        """distinct supports of tight forms (|D| = r >= 1)."""
        idx = np.nonzero((self.slack == 0) & (self.r >= 1))[0]
        return sorted(set(int(self.Dmask[q]) for q in idx))

    def zero_block(self):
        """rm-theory 8.7, my implementation.  Returns (excluded?, certificate or None, max clique size).
        joined[i][j] = some tight support contains i and j (i = j allowed: 'live').  For every form F with
        |D_F| + r_F < 2*|Y| for its largest clique Y inside D_F, the instance is excluded."""
        M = self.M; ts = self.tight_supports()
        J = [[False] * M for _ in range(M)]
        for d in ts:
            mem = [i for i in range(M) if d >> i & 1]
            for i in mem:
                for j in mem: J[i][j] = True
        live = [i for i in range(M) if J[i][i]]
        # largest clique overall (Theorem 1 check)
        best_all = max_clique(live, J)
        if len(best_all) > self.W:
            return True, {"kind": "C", "Y": best_all}, len(best_all)
        c = len(best_all)
        # only forms with |D_F| + r_F < 2c can be excluded
        cand = np.nonzero((self.cnt + self.r) < 2 * c)[0]
        for q in cand:
            d = int(self.Dmask[q]); mem = [i for i in live if d >> i & 1]
            Y = max_clique(mem, J)
            if int(self.cnt[q]) + int(self.r[q]) < 2 * len(Y):
                pairs = [[i, j, next(t for t in ts if t >> i & 1 and t >> j & 1)] for i in Y for j in Y if i <= j]
                return True, {"kind": "Z", "F": int(q), "DF": d, "rF": int(self.r[q]), "Y": Y, "pairs": pairs}, c
        return False, None, c


def max_clique(verts, J):
    best = []
    vs = list(verts)
    def rec(cur, rest):
        nonlocal best
        if len(cur) > len(best): best = cur[:]
        for k, v in enumerate(rest):
            if all(J[v][u] for u in cur):
                rec(cur + [v], rest[k + 1:])
    rec([], vs)
    return best
