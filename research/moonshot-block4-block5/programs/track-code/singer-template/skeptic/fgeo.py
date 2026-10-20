"""fgeo.py -- ms-fano: the Fano-gadget geometry in the ENGINE's own letters (dir z = (1+i)/2 I + (1-i)/2 X_z).
Positions: a term g (a product of letters) sits at the coset B'g (B' = monomial letters, which keep L0); its frame is
g^-1 L0.  Distance of a term from the raw input = log2 row support of g; distance to the kernel's coset = log2 row
support of K g^-1.  The interval = terms on geodesics (d0 + dK = m).  For an interval term g, A(g) = {z : g Z_z g^-1 is
diagonal} (= U cap L0 as z-labels); its annihilator s(g) is the subspace of F_2^m the term has 'raised' (pf-duality
Lemma I1).  Engine free letters used here: affine address maps, Clifford diagonal phases i^(a.x + 2 x^T B x), Paulis."""
import sys, os, os, itertools, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../common'))
import numpy as np
from bm import Game

def bits(x, m):
    return [(x >> j) & 1 for j in range(m)]

def dot(a, b):
    return bin(a & b).count('1') & 1

class Geo:
    def __init__(self, m, seed=1):
        self.m = m; self.G = G = Game(m); self.n = n = 1 << m
        self.rng = random.Random(seed)
        self.Z = {z: np.diag([(-1) ** dot(z, x) for x in range(n)]).astype(complex) for z in range(n)}

    # ---- engine free letters
    def affine(self, A, c):
        """address map x -> A x + c (A as list of column images of e_j)"""
        m, n = self.m, self.n
        sig = []
        for x in range(n):
            y = c
            for j in range(m):
                if (x >> j) & 1:
                    y ^= A[j]
            sig.append(y)
        assert len(set(sig)) == n
        return sig

    def rand_GL(self):
        m = self.m
        while True:
            A = [self.rng.randrange(1, self.n) for _ in range(m)]
            # independent?
            span = {0}
            ok = True
            for a in A:
                if a in span:
                    ok = False; break
                span |= {s ^ a for s in span}
            if ok:
                return A

    def clifford_phase(self, a, B):
        """q(x) = sum a_j x_j + 2 sum_{j<k} B_jk x_j x_k  (mod 4)"""
        m, n = self.m, self.n
        q = []
        for x in range(n):
            b = bits(x, m)
            v = sum(a[j] * b[j] for j in range(m))
            v += 2 * sum(B.get((j, k), 0) * b[j] * b[k] for j in range(m) for k in range(j + 1, m))
            q.append(v % 4)
        return q

    def rand_free(self, pauli=True):
        m, n = self.m, self.n
        A = self.rand_GL(); c = self.rng.randrange(n)
        sig = self.affine(A, c)
        a = [self.rng.randrange(4) for _ in range(m)]
        B = {(j, k): self.rng.randrange(2) for j in range(m) for k in range(j + 1, m)}
        q = self.clifford_phase(a, B)
        if pauli:
            z = self.rng.randrange(n); q = [(q[x] + 2 * dot(z, x)) % 4 for x in range(n)]
        return self.G.mono(sig, q)

    # ---- positions
    def d0(self, g):
        return int(round(np.log2(self.G.rowsupp(g).max())))

    def dK(self, g):
        h = self.G.mul(self.G.K, np.conj(g).T)
        return int(round(np.log2(self.G.rowsupp(h).max())))

    def A_of(self, g):
        """{z : g Z_z g^-1 diagonal}"""
        out = []
        gi = np.conj(g).T
        for z in range(self.n):
            M = g @ self.Z[z] @ gi
            if np.abs(M - np.diag(np.diag(M))).max() < 1e-9:
                out.append(z)
        return out

    def s_of(self, g):
        """annihilator of A(g): the raised subspace, as a sorted tuple of its nonzero vectors"""
        A = self.A_of(g)
        return tuple(v for v in range(1, self.n) if all(dot(v, z) == 0 for z in A))

    def fkey(self, g):
        return self.G.frame_key(g)

def span_of(vs):
    S = {0}
    for v in vs:
        S |= {s ^ v for s in S}
    return S

if __name__ == '__main__':
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    geo = Geo(m); G = geo.G
    D1 = G.D1
    # rank-1 interval terms D1 b over raw input
    found = {}
    tries = 0
    while len(found) < (1 << m) - 1 and tries < 20000:
        tries += 1
        b = geo.rand_free()
        g = G.mul(D1, b)
        if geo.d0(g) == 1 and geo.dK(g) == m - 1:
            s = geo.s_of(g)
            if s not in found:
                found[s] = b
    print(f"m={m}: rank-1 interval directions found: {sorted(found)} after {tries} random free letters")
    # all neighbours of L0 (for the count)
    nb = set()
    for _ in range(3000):
        b = geo.rand_free(); g = G.mul(D1, b)
        nb.add(geo.fkey(g))
    print(f"distinct frames D1 b over L0 (all neighbours): {len(nb)} (expected 2(2^m-1) = {2*((1<<m)-1)})")
