"""vfg.py (v-fp-build): my own frame-game arithmetic, independent of fpa.py.
Elements of Sp(2m,2) are tuples of COLUMN images (g[j] = image of basis vector e_j as a bitmask; bits 0..m-1 are
the x-part a, bits m..2m-1 the z-part b). omega((a|b),(a'|b')) = a.b' + b.a'. L0 = span of the z-basis vectors.
Position of a walker lam: lam^-1 L0 (canonical RREF). Level: distance of lam^-1 L0 from L0.
"""
from fractions import Fraction as Fr

def pc(x): return bin(x).count("1") & 1

class FG:
    def __init__(s, m):
        s.m = m; s.N = 2 * m
        s.I = tuple(1 << j for j in range(s.N))
        s.L0 = s.rref([1 << (m + j) for j in range(m)])
        s.K = s.I
        for j in range(m): s.K = s.mul(s.tv(1 << j), s.K)
    def app(s, g, v):
        out = 0; j = 0
        while v:
            if v & 1: out ^= g[j]
            v >>= 1; j += 1
        return out
    def mul(s, g, h): return tuple(s.app(g, h[j]) for j in range(s.N))
    def inv(s, g):
        # Gaussian elimination on [g | I] by columns: solve g x = e_j for each j
        N = s.N
        # rows of g: row i bitmask over columns
        rows = [sum(((g[j] >> i) & 1) << j for j in range(N)) | (1 << (N + i)) for i in range(N)]
        for col in range(N):
            piv = next(r for r in range(col, N) if (rows[r] >> col) & 1)
            rows[col], rows[piv] = rows[piv], rows[col]
            for r in range(N):
                if r != col and (rows[r] >> col) & 1: rows[r] ^= rows[col]
        # now rows[i] = e_i | (row i of g^-1) << N
        ginv_rows = [rows[i] >> N for i in range(N)]
        return tuple(sum(((ginv_rows[i] >> j) & 1) << i for i in range(N)) for j in range(N))
    def omega(s, x, y):
        m = s.m; M = (1 << m) - 1
        return pc((x & M) & (y >> m)) ^ pc((x >> m) & (y & M))
    def tv(s, u):
        """transvection x -> x + omega(x,u) u"""
        return tuple((1 << j) ^ (u if s.omega(1 << j, u) else 0) for j in range(s.N))
    def rref(s, vecs):
        basis = []
        for v in vecs:
            for b in basis:
                if v ^ b < v: v ^= b  # reduce by leading bit
            if v:
                basis.append(v); basis.sort(reverse=True)
        # full reduction
        basis.sort(reverse=True)
        for i in range(len(basis)):
            hb = basis[i].bit_length() - 1
            for k in range(len(basis)):
                if k != i and (basis[k] >> hb) & 1: basis[k] ^= basis[i]
        return tuple(sorted(basis, reverse=True))
    def rank(s, vecs): return len(s.rref(vecs))
    def img(s, g, sub): return s.rref([s.app(g, v) for v in sub])
    def pos(s, lam): return s.img(s.inv(lam), s.L0)
    def dist(s, X, Y): return s.rank(list(X) + list(Y)) - s.m
    def level(s, lam): return s.dist(s.L0, s.pos(lam))
    def inP(s, g): return s.img(g, s.L0) == s.L0
    def from_rows(s, R):
        """convert fpa's row-bitmask matrix into my column format"""
        N = s.N
        return tuple(sum(((R[i] >> j) & 1) << i for i in range(N)) for j in range(N))
    # contents: dict (input, element) -> coefficient
    def lmul(s, g, c):
        out = {}
        for (l, h), v in c.items():
            k = (l, s.mul(g, h)); out[k] = out.get(k, 0) + v
        return out
def cadd(*pairs):
    out = {}
    for co, c in pairs:
        for k, v in c.items(): out[k] = out.get(k, 0) + co * v
    return out
