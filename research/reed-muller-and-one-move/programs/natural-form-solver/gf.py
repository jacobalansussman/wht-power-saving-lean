"""gf.py: small GF(2^m) toolkit for rm-theory (field model of the Singer no-phase game).

Elements are ints 0..q-1 in the polynomial basis of a primitive polynomial. Everything is exact (bit arithmetic).
Coordinates of the game: a SELF-DUAL basis f_1..f_m of F_q over F_2 (Tr(f_i f_j) = delta_ij), so that the
standard dot product of coordinate vectors is the trace form Tr(xy), and the target I is the trace form.
"""
import itertools, random

PRIM = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011, 8: 0b100011101,
        9: 0b1000010001, 10: 0b10000001001}


class GF:
    def __init__(self, m, poly=None):
        self.m = m; self.q = 1 << m; self.n = self.q - 1
        self.poly = poly if poly is not None else PRIM[m]
        # log / exp tables with alpha = x (primitive)
        self.exp = [0] * (2 * self.n); self.log = [None] * self.q
        x = 1
        for k in range(self.n):
            self.exp[k] = x; self.log[x] = k
            x <<= 1
            if x & self.q: x ^= self.poly
        assert x == 1 and len(set(self.exp[:self.n])) == self.n, "polynomial not primitive"
        for k in range(self.n, 2 * self.n): self.exp[k] = self.exp[k - self.n]
        self.trv = [self._tr(a) for a in range(self.q)]

    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        return self.exp[self.log[a] + self.log[b]]

    def pw(self, a, e):
        if a == 0: return 0 if e > 0 else 1
        return self.exp[(self.log[a] * e) % self.n]

    def inv(self, a):
        return self.exp[(-self.log[a]) % self.n]

    def _tr(self, a):
        s = 0; x = a
        for _ in range(self.m):
            s ^= x; x = self.mul(x, x)
        assert s in (0, 1)
        return s

    def tr(self, a):
        return self.trv[a]

    def alpha(self, k):
        return self.exp[k % self.n]

    def self_dual_basis(self, seed=1):
        """Random greedy search for f_1..f_m with Tr(f_i f_j) = delta_ij."""
        rng = random.Random(seed)
        elts = list(range(1, self.q))
        for attempt in range(100000):
            basis = []
            ok = True
            for j in range(self.m):
                cands = [x for x in elts if self.tr(self.mul(x, x)) == 1 and all(self.tr(self.mul(x, f)) == 0 for f in basis)]
                # x must also be independent of basis: orthonormality forces it
                if not cands: ok = False; break
                basis.append(rng.choice(cands))
            if ok: return basis
        raise RuntimeError("no self-dual basis found")

    def coords(self, x, basis):
        """coordinates in a self-dual basis: x_j = Tr(x f_j)."""
        return tuple(self.tr(self.mul(x, f)) for f in basis)


def gold_exponents(m):
    """exponents 1 + 2^i, i = 1..floor(m/2): the cyclotomic classes of RM(2,m) quadratic part."""
    return [1 + (1 << i) for i in range(1, m // 2 + 1)]
