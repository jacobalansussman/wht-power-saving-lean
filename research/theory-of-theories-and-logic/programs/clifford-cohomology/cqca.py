# cqca.py -- exact Pauli / Clifford / Clifford-QCA arithmetic (pure python, no imports beyond random).
# meta-cohomology, 14 Oct 2026.  Conventions (all exact, phases mod 4):
#   A Pauli on n qubits is (x, z, k) = i^k X^x Z^z, x and z bit masks (python ints, unbounded).
#   A Clifford map c is the conjugation action P -> U P U^dagger; it is stored by the images of the generators.
#   One-column Clifford (m qubits): list of 2m Paulis, images of X_0..X_{m-1}, Z_0..Z_{m-1}.
#   Clifford QCA (bulk rule) on the chain Z, m qubits per site: list of 2m Paulis on the chain, images of the
#   generators at site 0; qubit (site s, row q) has bit index (s + OFF) * m + q.  The image of a generator at
#   site s is the image at site 0 shifted by s columns (translation invariance).
import random

OFF = 400          # site offset for chain Paulis; all supports stay far inside [-OFF, OFF]


def pc(a):
    return bin(a).count("1")


def pmul(P, Q):
    """(i^a X^x1 Z^z1)(i^b X^x2 Z^z2) = i^(a+b) (-1)^|z1 & x2| X^(x1^x2) Z^(z1^z2)."""
    x1, z1, a = P
    x2, z2, b = Q
    return (x1 ^ x2, z1 ^ z2, (a + b + 2 * pc(z1 & x2)) % 4)


def herm(x, z):
    """Hermitian Pauli with labels (x, z): i^|x&z| X^x Z^z."""
    return (x, z, pc(x & z) % 4)


def is_herm(P):
    x, z, k = P
    return (k - pc(x & z)) % 2 == 0


def symp(P, Q):
    """symplectic form <P, Q> in F_2 (1 iff they anticommute)."""
    return (pc(P[0] & Q[1]) + pc(P[1] & Q[0])) % 2


# ---------------------------------------------------------------- Clifford maps on a fixed set of n qubits
def capply(imgX, imgZ, P):
    """Apply a Clifford map given by images of X_b (dict or list) and Z_b to a Pauli P (phases exact)."""
    x, z, k = P
    R = (0, 0, k)
    b = 0
    while x:
        if x & 1:
            R = pmul(R, imgX(b))
        x >>= 1
        b += 1
    b = 0
    while z:
        if z & 1:
            R = pmul(R, imgZ(b))
        z >>= 1
        b += 1
    return R


class Cliff:
    """Clifford map on n qubits (n = len(X)), images of X_b, Z_b."""
    def __init__(self, X, Z):
        self.X = list(X)
        self.Z = list(Z)
        self.n = len(X)

    @staticmethod
    def ident(n):
        return Cliff([(1 << b, 0, 0) for b in range(n)], [(0, 1 << b, 0) for b in range(n)])

    def __call__(self, P):
        return capply(lambda b: self.X[b], lambda b: self.Z[b], P)

    def __mul__(self, other):          # (self * other)(P) = self(other(P))
        return Cliff([self(p) for p in other.X], [self(p) for p in other.Z])

    def __eq__(self, other):
        return self.X == other.X and self.Z == other.Z

    def key(self):
        return (tuple(self.X), tuple(self.Z))

    def valid(self):
        """images Hermitian and with the canonical commutation relations."""
        G = self.X + self.Z
        n = self.n
        for i, P in enumerate(G):
            if not is_herm(P):
                return False
            for j, Q in enumerate(G):
                want = 1 if (abs(i - j) == n) else 0
                if symp(P, Q) != want:
                    return False
        return True

    def is_identity(self):
        return self == Cliff.ident(self.n)

    def label_matrix_cols(self):
        """columns of the label matrix: for generator g (X_0..X_{n-1}, Z_0..Z_{n-1}) the label (x, z) of its image."""
        return [(p[0], p[1]) for p in self.X + self.Z]


# ---------------------------------------------------------------- gates as conjugation rules on single Paulis
def g_H(q):
    def f(P):
        x, z, k = P
        bx, bz = (x >> q) & 1, (z >> q) & 1
        # H X H = Z, H Z H = X, H Y H = -Y.  Write P = i^k X^x Z^z; factor on qubit q is X^bx Z^bz.
        # conj(X^bx Z^bz) = Z^bx X^bz = (-1)^(bx bz) X^bz Z^bx
        x2 = (x & ~(1 << q)) | (bz << q)
        z2 = (z & ~(1 << q)) | (bx << q)
        # careful: the factor on qubit q sits inside X^x Z^z; reordering X^x Z^z -> ... only qubit q changes,
        # other qubits' factors commute with qubit q's factors.
        return (x2, z2, (k + 2 * (bx & bz)) % 4)
    return f


def g_S(q):
    def f(P):
        x, z, k = P
        bx = (x >> q) & 1
        # S X S^dag = Y = i X Z ; S Z S^dag = Z.  X^bx Z^bz -> (i X Z)^bx Z^bz = i^bx X^bx Z^(bx+bz)
        if not bx:
            return P
        # moving the new Z_q into Z^z: X^x Z^z with an extra Z_q placed right after X_q: Z_q commutes with
        # X on other qubits, so X^x Z_q Z^z = X^x Z^(z ^ e_q) (Z's commute)
        return (x, z ^ (1 << q), (k + 1) % 4)
    return f


def g_CNOT(c, t):
    def f(P):
        x, z, k = P
        bxc = (x >> c) & 1
        bzt = (z >> t) & 1
        # X_c -> X_c X_t ; Z_t -> Z_c Z_t ; X_t, Z_c fixed.  On X^x Z^z: X part gains X_t^bxc, Z part gains Z_c^bzt.
        # X^x' Z^z' with x' = x ^ (bxc << t), z' = z ^ (bzt << c); no reordering needed (X's commute, Z's commute).
        return (x ^ (bxc << t), z ^ (bzt << c), k)
    return f


def g_CZ(a, b):
    def f(P):
        x, z, k = P
        bxa = (x >> a) & 1
        bxb = (x >> b) & 1
        # X_a -> X_a Z_b, X_b -> Z_a X_b.  conj(X^x Z^z) = (prod over X's of images) Z^z.
        # images: X_a Z_b and Z_a X_b; product of X-images in ascending bit order then Z^z.
        R = (0, 0, k)
        xx = x
        bit = 0
        while xx:
            if xx & 1:
                if bit == a:
                    R = pmul(R, (1 << a, 1 << b, 0))
                elif bit == b:
                    # Z_a X_b with a != b: they commute, so Z_a X_b = X_b Z_a = (x=e_b, z=e_a, phase 0)
                    R = pmul(R, (1 << b, 1 << a, 0))
                else:
                    R = pmul(R, (1 << bit, 0, 0))
            xx >>= 1
            bit += 1
        return pmul(R, (0, z, 0))
    return f


def g_X(q):
    def f(P):
        x, z, k = P
        return (x, z, (k + 2 * ((z >> q) & 1)) % 4)
    return f


def g_Z(q):
    def f(P):
        x, z, k = P
        return (x, z, (k + 2 * ((x >> q) & 1)) % 4)
    return f


def g_rot(Qx, Qz):
    """conjugation by D = (1+i)/2 + (1-i)/2 Q, Q the Hermitian Pauli with labels (Qx, Qz):
    P -> P if P commutes with Q, else P -> i P Q.  (D = sqrt(Q) up to a phase: a direction D_z when Q = X_z.)"""
    Q = herm(Qx, Qz)

    def f(P):
        if (pc(P[0] & Qz) + pc(P[1] & Qx)) % 2 == 0:
            return P
        return pmul((0, 0, 1), pmul(P, Q))
    return f


def apply_gates(gates, P):
    """conjugate P by the unitary U = g_last ... g_first (gates listed in time order)."""
    for g in gates:
        P = g(P)
    return P


# ---------------------------------------------------------------- chain (QCA) helpers
def bit(m, s, q):
    return (s + OFF) * m + q


def shiftP(P, s, m):
    x, z, k = P
    if s >= 0:
        return (x << (s * m), z << (s * m), k)
    return (x >> (-s * m), z >> (-s * m), k)


def sites_of(P, m):
    x, z, _ = P
    a = x | z
    out = []
    b = 0
    while a:
        if a & 1:
            out.append(b // m - OFF)
        a >>= 1
        b += 1
    return sorted(set(out))


def site_factor(P, m, s):
    """labels (x_s, z_s) of P on site s, as m-bit masks."""
    sh = (s + OFF) * m
    mask = (1 << m) - 1
    return ((P[0] >> sh) & mask, (P[1] >> sh) & mask)


class QCA:
    """Clifford QCA on the chain: rule[i] = image of generator i at site 0 (i < m: X_i; i >= m: Z_{i-m})."""
    def __init__(self, m, rule):
        self.m = m
        self.rule = list(rule)

    @staticmethod
    def ident(m):
        return QCA(m, [(1 << bit(m, 0, q), 0, 0) for q in range(m)] + [(0, 1 << bit(m, 0, q), 0) for q in range(m)])

    def img_gen(self, b, isz):
        s, q = b // self.m - OFF, b % self.m
        return shiftP(self.rule[self.m * isz + q], s, self.m)

    def __call__(self, P):
        return capply(lambda b: self.img_gen(b, 0), lambda b: self.img_gen(b, 1), P)

    def __mul__(self, other):          # (self * other)(P) = self(other(P))
        return QCA(self.m, [self(p) for p in other.rule])

    def __eq__(self, other):
        return self.rule == other.rule

    def radius(self):
        r = 0
        for p in self.rule:
            ss = sites_of(p, self.m)
            if ss:
                r = max(r, max(abs(s) for s in ss))
        return r

    def symbol(self):
        """dict d -> list of 2m columns, column i = label (x_d, z_d) (m-bit masks) of the site-d factor of image i."""
        out = {}
        for i, p in enumerate(self.rule):
            for s in sites_of(p, self.m):
                out.setdefault(s, [(0, 0)] * (2 * self.m))
                col = list(out[s])
                col[i] = site_factor(p, self.m, s)
                out[s] = col
        return out

    def symbol_is_identity(self):
        S = self.symbol()
        if set(S.keys()) - {0}:
            return False
        I = [(1 << q, 0) for q in range(self.m)] + [(0, 1 << q) for q in range(self.m)]
        return S.get(0, I) == I


def layer_qca(m, gate_makers, reach=6):
    """QCA of a translation-invariant layer: gate_makers(j) -> list of gate functions for the gates attached to
    site j (all gates of the layer must commute).  Image of a generator = conjugation by all gates near it."""
    rule = []
    for isz in (0, 1):
        for q in range(m):
            P = (1 << bit(m, 0, q), 0, 0) if isz == 0 else (0, 1 << bit(m, 0, q), 0)
            for j in range(-reach, reach + 1):
                for g in gate_makers(j):
                    P = g(P)
            rule.append(P)
    return QCA(m, rule)


# ---------------------------------------------------------------- the two candidate splittings
def collapse(Q):
    """Pi_0: replace the tensor product over sites by the product over sites in ascending order, on one column."""
    m = Q.m
    imgs = []
    for p in Q.rule:
        R = (0, 0, p[2])
        for s in sites_of(p, m):
            xs, zs = site_factor(p, m, s)
            R = pmul(R, (xs, zs, 0))
        imgs.append(R)
    return Cliff(imgs[:m], imgs[m:])


def ring_restrict(Q, L):
    """Psi_L: wrap the rule on a ring of L sites (L odd, L > 2*radius+1); apply it to the logical Paulis
    X_q^(x)L, Z_q^(x)L (product of the images of X_{q,j}, j = 0..L-1, in that order); read the result as
    i^k X^(a..a) Z^(b..b) and return the one-column image (a, b, k).  Returns None if not translation invariant."""
    m = Q.m
    assert L % 2 == 1

    def wrap(P):
        x, z, k = P
        xr = zr = 0
        for s in sites_of(P, m):
            xs, zs = site_factor(P, m, s)
            t = s % L
            assert ((xr >> (t * m)) | (zr >> (t * m))) & ((1 << m) - 1) == 0, "collision: L too small"
            xr |= xs << (t * m)
            zr |= zs << (t * m)
        return (xr, zr, k)

    imgs = []
    for i, p in enumerate(Q.rule):
        R = (0, 0, 0)
        for j in range(L):
            R = pmul(R, wrap(shiftP(p, j, m)))
        x, z, k = R
        mask = (1 << m) - 1
        a, b = x & mask, z & mask
        rep_a = sum(a << (t * m) for t in range(L))
        rep_b = sum(b << (t * m) for t in range(L))
        if x != rep_a or z != rep_b:
            return None
        imgs.append((a, b, k))
    return Cliff(imgs[:m], imgs[m:])


def gamma_matrix(Q):
    """gamma(e_i, e_j) = sum_{d > d'} <M_d e_i, M_d' e_j> for the 2m label basis vectors; returns the 2m x 2m F2 matrix."""
    S = Q.symbol()
    ds = sorted(S.keys())
    n = 2 * Q.m
    G = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            t = 0
            for d in ds:
                for d2 in ds:
                    if d > d2:
                        a = S[d][i]
                        b = S[d2][j]
                        t ^= (pc(a[0] & b[1]) + pc(a[1] & b[0])) & 1
            G[i][j] = t
    return G


def chi_form(Q):
    """chi(e_i) = <M_odd e_i, M_even e_i> (the predicted difference between the ring splitting and the collapse)."""
    S = Q.symbol()
    n = 2 * Q.m
    out = []
    for i in range(n):
        ox = oz = ex = ez = 0
        for d, cols in S.items():
            x, z = cols[i]
            if d % 2:
                ox ^= x; oz ^= z
            else:
                ex ^= x; ez ^= z
        out.append((pc(ox & ez) + pc(oz & ex)) & 1)
    return out
