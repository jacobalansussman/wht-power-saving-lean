"""vfield.py (v-rm2-solve): independent GF(2^m) model. Targets as symmetric F2 matrices (bitmask of upper triangle),
directions as ints 1..2^m-1 in a self-dual basis found by exhaustive search (not greedy-random)."""
import itertools
PRIM = {3: 0b1011, 4: 0b10011, 5: 0b100101}
class Fq:
    def __init__(s, m):
        s.m = m; s.q = 1 << m; s.n = s.q - 1; p = PRIM[m]
        s.e = []; x = 1
        for _ in range(s.n):
            s.e.append(x); x <<= 1
            if x & s.q: x ^= p
        assert x == 1 and len(set(s.e)) == s.n
        s.lg = {v: i for i, v in enumerate(s.e)}
    def mul(s, a, b):
        if a == 0 or b == 0: return 0
        return s.e[(s.lg[a] + s.lg[b]) % s.n]
    def tr(s, a):
        t = 0; x = a
        for _ in range(s.m): t ^= x; x = s.mul(x, x)
        assert t in (0, 1); return t
    def sdbasis(s):
        # lexicographically first self-dual basis (exhaustive, deterministic)
        els = range(1, s.q)
        def rec(B):
            if len(B) == s.m: return B
            for x in els:
                if (not B or x > B[-1]) and s.tr(s.mul(x, x)) == 1 and all(s.tr(s.mul(x, f)) == 0 for f in B):
                    r = rec(B + [x])
                    if r: return r
            return None
        return rec([])
class Model:
    def __init__(s, m):
        s.F = F = Fq(m); s.m = m; s.B = F.sdbasis()
        s.mons = [(i, j) for i in range(m) for j in range(i, m)]
    def mulmat(s, a):
        F = s.F; m = s.m
        return [[F.tr(F.mul(F.mul(s.B[i], a), s.B[j])) for j in range(m)] for i in range(m)]
    def target(s, k):  # T_k = multiplication by alpha^(2k)
        return s.mulmat(s.F.e[(2 * k) % s.F.n])
    def smask(s, T):
        return sum((T[i][j] & 1) << q for q, (i, j) in enumerate(s.mons))
    def dmask(s, z):
        b = [(z >> j) & 1 for j in range(s.m)]
        return sum((b[i] & b[j]) << q for q, (i, j) in enumerate(s.mons))
    def Qval(s, c, z):
        # form with coefficient bitmask c over monomials; Q(z) = <c, dmask(z)>
        return bin(c & s.dmask(z)).count("1") & 1
    def par(s, c, T):
        return bin(c & s.smask(T)).count("1") & 1
