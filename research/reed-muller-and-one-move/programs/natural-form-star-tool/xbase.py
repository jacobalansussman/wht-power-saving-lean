"""xbase.py (rm2-exact, 15 Oct 2026): own field / forms / targets code for the multi-target system (Theorem B form).

Coordinates: a self-dual basis f_1..f_m of GF(2^m) (poly x^5+x^2+1 at m=5, x^4+x+1 at m=4, x^3+x+1 at m=3), so the
dot product is Tr(xy).  Directions z = nonzero coordinate vectors (ints 1..2^m-1, bit j = coordinate j).
Forms: Q(z) = sum_{j<=k} q_jk z_j z_k, q packed as an int over the pairs (j,k), j<=k (index order PAIRS).
Pairing with a symmetric matrix S: <Q,S> = sum_{j<=k} q_jk S_jk ; Q(z) = <Q, z z^T>.
Targets: canonical Singer targets T_k = mult(alpha^(2k)) (copy k of rm-theory), or any symmetric matrix.
Everything exact (ints / bit masks).
"""
import itertools
import numpy as np

PRIM = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011}


class Field:
    def __init__(self, m):
        self.m = m; self.q = 1 << m; self.n = self.q - 1
        poly = PRIM[m]
        self.exp = [0] * (2 * self.n); self.log = [None] * self.q
        x = 1
        for k in range(self.n):
            self.exp[k] = x; self.log[x] = k
            x <<= 1
            if x & self.q: x ^= poly
        assert x == 1
        for k in range(self.n, 2 * self.n): self.exp[k] = self.exp[k - self.n]
        self.trv = []
        for a in range(self.q):
            s, y = 0, a
            for _ in range(m): s ^= y; y = self.mul(y, y)
            assert s in (0, 1); self.trv.append(s)
        self.basis = self._self_dual()

    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        return self.exp[self.log[a] + self.log[b]]

    def tr(self, a): return self.trv[a]

    def _self_dual(self):
        # deterministic backtracking for an orthonormal basis of the trace form
        m = self.m
        cand = [x for x in range(1, self.q) if self.tr(self.mul(x, x)) == 1]
        def bt(B):
            if len(B) == m: return B
            for x in cand:
                if (not B or x > B[-1]) and all(self.tr(self.mul(x, f)) == 0 for f in B):
                    r = bt(B + [x])
                    if r: return r
            return None
        B = bt([]); assert B
        return B

    def coord(self, x):
        return sum(self.tr(self.mul(x, f)) << j for j, f in enumerate(self.basis))

    def mult_matrix(self, g):
        m = self.m
        return [[self.tr(self.mul(g, self.mul(self.basis[j], self.basis[k]))) for k in range(m)] for j in range(m)]


class Forms:
    def __init__(self, m):
        self.m = m
        self.PAIRS = [(j, k) for j in range(m) for k in range(j, m)]
        self.N = len(self.PAIRS); self.nf = 1 << self.N
        self.dirs = list(range(1, 1 << m))            # move i <-> direction dirs[i]
        nd = len(self.dirs)
        mon = np.zeros((nd, self.N), dtype=np.int64)
        for i, z in enumerate(self.dirs):
            for p, (j, k) in enumerate(self.PAIRS):
                mon[i, p] = ((z >> j) & 1) & ((z >> k) & 1)
        self.mon = mon
        qs = np.arange(self.nf, dtype=np.int64)
        qbits = ((qs[:, None] >> np.arange(self.N)[None, :]) & 1)          # nf x N
        vals = (qbits @ mon.T) & 1                                          # nf x nd : Q(z_i)
        self.vals = vals.astype(np.uint8)
        self.mask = (vals.astype(np.int64) << np.arange(nd)[None, :]).sum(axis=1)   # support bit masks
        self.wt = vals.sum(axis=1)
        self.qbits = qbits

    def tvec(self, S):
        """upper-triangle vector (int) of a symmetric 0/1 matrix S"""
        return sum(int(S[j][k]) << p for p, (j, k) in enumerate(self.PAIRS))

    def parity(self, tv):
        """array over all forms of <Q, T> for target vector tv"""
        tb = np.array([(tv >> p) & 1 for p in range(self.N)], dtype=np.int64)
        return (self.qbits @ tb) & 1


def canon_targets(Fd, ks):
    """canonical Singer targets T_k = mult(alpha^(2k)) for copy indices ks"""
    return [Fd.mult_matrix(Fd.exp[(2 * k) % Fd.n]) for k in ks]


def identity(m): return [[int(j == k) for k in range(m)] for j in range(m)]


class System:
    """one move per listed direction multiset (counts c[z]), W arrays with targets; per-form r, s."""
    def __init__(self, fo, targets, counts=None):
        self.fo = fo; self.W = len(targets)
        self.counts = counts if counts is not None else {z: 1 for z in fo.dirs}
        P = np.stack([fo.parity(fo.tvec(T)) for T in targets], axis=1)   # nf x W
        self.P = P
        self.r = P.sum(axis=1)
        cvec = np.array([self.counts.get(z, 0) for z in fo.dirs], dtype=np.int64)
        self.cvec = cvec
        self.d = fo.vals.astype(np.int64) @ cvec          # |D_Q| (moves)
        self.s = self.d - self.r                          # rank E[D,D] exactly
    def violations(self):
        return np.nonzero(self.s < 0)[0]
    def tight(self):
        return np.nonzero((self.s == 0) & (self.r > 0))[0]


if __name__ == "__main__":
    for m in (3, 4, 5):
        Fd = Field(m); fo = Forms(m)
        # checks: sigma symmetric, mult matrices symmetric, trace parity tr Q = <Q, I>
        a = Fd.exp[1]; S = Fd.mult_matrix(a)
        assert all(S[j][k] == S[k][j] for j in range(m) for k in range(m))
        I = identity(m); pI = fo.parity(fo.tvec(I))
        trq = sum(fo.vals[:, fo.dirs.index(1 << j)].astype(int) for j in range(m)) & 1
        assert (pI == trq).all()
        # T_0 = mult(1) = I
        assert canon_targets(Fd, [0])[0] == I
        # canonical: copy k odd at Q exactly when ... n_odd(Q) in {0, 2^(m-1)} over all n copies (rm-theory F3)
        allT = canon_targets(Fd, range(Fd.n))
        sy = System(fo, allT)
        vals = set(sy.r.tolist())
        print("m=%d self-dual basis %s; n_odd values over all %d copies: %s; forms %d" % (m, Fd.basis, Fd.n, sorted(vals), fo.nf))
