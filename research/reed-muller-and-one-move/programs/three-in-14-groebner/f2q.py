"""f2q.py (r3-3in14): my own small toolkit, written from the definitions (no code shared with rm2-solve or
v-rm2-solve).

- GF(2^m) by log/exp tables, alpha a root of the same primitive polynomials (x^3+x+1, x^4+x+1, x^5+x^2+1).
- A self-dual basis f_1..f_m (Tr(f_i f_j) = delta_ij), found by a deterministic depth-first search (the FIRST one in
  increasing integer order, so in general a different basis from the other tracks' random one).
- Coordinates of x: (Tr(x f_1), ..., Tr(x f_m)), packed as an int (bit j = coordinate j).
- Canonical Singer target T_k: the symmetric matrix (Tr(alpha^(2k) f_i f_j))_ij  (= sum over the frame alpha^k f_j).
- Binary quadratic forms Q(z) = sum_{i<=j} q_ij z_i z_j, as bitmasks over the monomials (i, j), i <= j.
  Q(z) = <Q, z z^T> and p_T(Q) = <Q, T> = sum_{i<=j} q_ij T_ij (mod 2).
"""
import itertools

PRIM = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101}


def popc(x):
    return bin(x).count("1")


class Field:
    def __init__(self, m):
        self.m = m; self.q = 1 << m; self.n = self.q - 1; poly = PRIM[m]
        self.exp = [0] * (2 * self.n); self.log = {}
        x = 1
        for k in range(self.n):
            self.exp[k] = x; self.log[x] = k
            x <<= 1
            if x & self.q: x ^= poly
        assert x == 1 and len(self.log) == self.n
        for k in range(self.n, 2 * self.n): self.exp[k] = self.exp[k - self.n]
        self.basis = self._self_dual()
        self.mons = [(i, i) for i in range(m)] + [(i, j) for i in range(m) for j in range(i + 1, m)]
        self.P = len(self.mons)
        self.coord = {x: sum(self.tr(self.mul(x, f)) << j for j, f in enumerate(self.basis)) for x in range(self.q)}
        assert len(set(self.coord.values())) == self.q        # coordinates are a bijection

    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        return self.exp[self.log[a] + self.log[b]]

    def alpha(self, k):
        return self.exp[k % self.n]

    def tr(self, a):
        s, x = 0, a
        for _ in range(self.m):
            s ^= x; x = self.mul(x, x)
        assert s in (0, 1)
        return s

    def _self_dual(self):
        m = self.m
        cand = [x for x in range(1, self.q) if self.tr(self.mul(x, x)) == 1]
        def rec(b):
            if len(b) == m: return b
            for x in cand:
                if (not b or x > b[-1]) and all(self.tr(self.mul(x, f)) == 0 for f in b):
                    r = rec(b + [x])
                    if r: return r
            return None
        b = rec([])
        assert b is not None
        return b

    # ---- symmetric matrices as bitmasks over monomials
    def sym_mask(self, S):
        return sum((S[i][j] & 1) << t for t, (i, j) in enumerate(self.mons))

    def dir_mask(self, z):
        b = [(z >> j) & 1 for j in range(self.m)]
        return sum((b[i] & b[j]) << t for t, (i, j) in enumerate(self.mons))

    def target(self, k):
        m = self.m; c = self.alpha(2 * k)
        return [[self.tr(self.mul(c, self.mul(self.basis[i], self.basis[j]))) for j in range(m)] for i in range(m)]

    def frame_canonical(self, k):
        return [self.coord[self.mul(self.alpha(k), f)] for f in self.basis]


def check_field(F):
    """basic identities: T_k symmetric, T_k = T_1^k over F2, canonical frame sums to T_k."""
    m = F.m
    def matmul(A, B):
        return [[sum(A[i][t] & B[t][j] for t in range(m)) & 1 for j in range(m)] for i in range(m)]
    T1 = F.target(1); P = [[int(i == j) for j in range(m)] for i in range(m)]
    for k in range(F.n):
        T = F.target(k)
        assert all(T[i][j] == T[j][i] for i in range(m) for j in range(m))
        assert T == P, (k,)
        P = matmul(P, T1)
        v = 0
        for z in F.frame_canonical(k): v ^= F.dir_mask(z)
        assert v == F.sym_mask(T)
    assert F.target(0) == [[int(i == j) for j in range(m)] for i in range(m)]
    return True


def all_frames(F, k):
    """all m-sets of distinct nonzero directions with sum z z^T = T_k."""
    goal = F.sym_mask(F.target(k)); dm = {z: F.dir_mask(z) for z in range(1, F.q)}
    out = []
    for S in itertools.combinations(range(1, F.q), F.m):
        v = 0
        for z in S: v ^= dm[z]
        if v == goal: out.append(S)
    return out


class Config:
    """W arrays with symmetric targets Ts (lists), moves = list of direction ints (with repeats; order irrelevant for
    every test here). For every form Q (all 2^P bitmasks): Dmask[Q] = bitmask of moves i with Q(z_i) = 1,
    pmask[Q] = bitmask of arrays l odd at Q."""
    def __init__(self, F, Ts, moves):
        self.F, self.Ts, self.moves = F, Ts, list(moves)
        self.W, self.M = len(Ts), len(moves)
        zm = [F.dir_mask(z) for z in self.moves]; tm = [F.sym_mask(T) for T in Ts]
        NQ = 1 << F.P
        self.Dmask = [0] * NQ; self.pmask = [0] * NQ
        for Q in range(NQ):
            d = 0
            for i, x in enumerate(zm):
                if popc(Q & x) & 1: d |= 1 << i
            p = 0
            for l, x in enumerate(tm):
                if popc(Q & x) & 1: p |= 1 << l
            self.Dmask[Q] = d; self.pmask[Q] = p

    def span_deletion(self):
        """(span_ok, deletion_ok): span = no form with empty support and an odd array; deletion = |D| >= r always."""
        span = all(not (self.Dmask[Q] == 0 and self.pmask[Q]) for Q in range(len(self.Dmask)))
        dele = all(popc(self.Dmask[Q]) >= popc(self.pmask[Q]) for Q in range(len(self.Dmask)))
        return span, dele

    def tight_supports(self):
        """set of supports D (bitmask) of tight forms: |D| = r >= 1. Returns dict D -> (Q, pmask)."""
        out = {}
        for Q in range(len(self.Dmask)):
            d, p = self.Dmask[Q], self.pmask[Q]
            if p and popc(d) == popc(p):
                if d in out:
                    assert out[d][1] == p, "same support, different parities: span fails"
                else:
                    out[d] = (Q, p)
        return out
