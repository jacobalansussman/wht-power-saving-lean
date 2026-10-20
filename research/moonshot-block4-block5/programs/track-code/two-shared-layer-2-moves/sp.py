"""sp.py -- w5-twoshared: frame-level arithmetic in the ENGINE game through the symplectic picture.
A Clifford element g acts on Paulis (x|z) in F2^8 (bits 0-3 = x, 4-7 = z) by S_g(P) = g P g^-1.  Its frame B'g
corresponds to the Lagrangian Lag(g) = S_g^-1(L0), L0 = {(0|z)}; its layer lambda(g) = 4 - dim(Lag(g) cap L0)
(row support 2^lambda); the kernel frame is Lag(K) = {(z|z)}.  Engine letters = Stab(L0) = {(x;z) -> (Px; P^-T z + B P x)},
B symmetric; D = dir(e1): x_1 += z_1.  S_(gh) = S_g S_h.  Maps are tuples of 8 images of the basis vectors."""
import random
def app(S, v):
    r = 0; i = 0
    while v:
        if v & 1: r ^= S[i]
        v >>= 1; i += 1
    return r
def comp(S, T):              # S o T
    return tuple(app(S, t) for t in T)
I8 = tuple(1 << i for i in range(8))
def omega(u, v):
    return (bin((u & 15) & (v >> 4)).count('1') + bin((u >> 4) & (v & 15)).count('1')) & 1
def inv(S):
    """symplectic inverse: S^-1(w) = sum_i omega(S(J e_i), w)... use: S^-1 = Omega S^T Omega"""
    # S^-1(w)_i : coefficient c_i with w = S(sum c_i e_i); omega(S e_j, w) = omega(e_j, S^-1 w)
    out = []
    for i in range(8):
        # image of e_i under S^-1: vector u with omega(e_j, u) = omega(S e_j, e_i) for all j
        u = 0
        for j in range(8):
            if omega(S[j], 1 << i):
                # omega(e_j, u) = 1: e_j x-bit j<4 pairs with z-bit j+4, and vice versa
                u |= 1 << ((j + 4) % 8)
        out.append(u)
    return tuple(out)
D = tuple((1 << 4) | 1 if i == 4 else (1 << i) for i in range(8))   # e_z1 -> e_z1 + e_x1
def span(vs):
    s = {0}
    for v in vs:
        if v not in s:
            s |= {w ^ v for w in s}
    return frozenset(s)
L0 = span([1 << 4, 1 << 5, 1 << 6, 1 << 7])
LK = span([(1 << 4) | 1, (1 << 5) | 2, (1 << 6) | 4, (1 << 7) | 8])
def lag(S):
    Si = inv(S)
    return span([app(Si, 1 << (4 + j)) for j in range(4)])
def lamL(L):
    return 4 - (len(L & L0).bit_length() - 1)
def lam(S):
    return lamL(lag(S))
def gl4(rng):
    while True:
        cols = [rng.randrange(1, 16) for _ in range(4)]
        if len(span(cols)) == 16: return cols
def letter(rng):
    """random engine letter: (x;z) -> (Px; P^-T z + B P x)"""
    P = gl4(rng)                      # P e_j = cols[j]
    # P^-T: for z, need Q with Q = (P^-1)^T: omega-dual; compute P^-1 by brute force on 16 vectors
    Pm = {}
    for v in range(16):
        w = 0
        for j in range(4):
            if (v >> j) & 1: w ^= P[j]
        Pm[v] = w
    Pinv = {w: v for v, w in Pm.items()}
    # (P^-T z)_i = sum_j (P^-1)_{j i} z_j  ->  image of e_zi under P^-T is sum_j (P^-1)_{j,i} e_zj, (P^-1)_{j,i} = bit j of P^-1 e_i
    PT = [0] * 4
    for i in range(4):
        col = Pinv[1 << i]            # P^-1 e_i, bits j = (P^-1)_{j,i}
        # P^-T e_i = sum_j (P^-T)_{j,i} e_j = sum_j (P^-1)_{i,j} e_j ; (P^-1)_{i,j} = bit i of P^-1 e_j
        PT[i] = sum(((Pinv[1 << j] >> i) & 1) << j for j in range(4))
    Bm = [[0] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(i, 4):
            Bm[i][j] = Bm[j][i] = rng.randrange(2)
    def Bx(x):
        r = 0
        for i in range(4):
            if sum(Bm[i][j] * ((x >> j) & 1) for j in range(4)) & 1: r |= 1 << i
        return r
    S = []
    for i in range(4):                # x basis e_i -> (P e_i ; B P e_i)
        px = P[i]; S.append(px | (Bx(px) << 4))
    for i in range(4):                # z basis -> (0 ; P^-T e_i)
        S.append(PT[i] << 4)
    S = tuple(S)
    assert lag(S) == L0
    return S
def is_symp(S):
    return all(omega(S[i], S[j]) == omega(1 << i, 1 << j) for i in range(8) for j in range(8))
def mv(a, g):
    """S of D a g"""
    return comp(D, comp(a, g))
def prog_letter(rng, gs, up=1, tries=2000):
    """letter a with lam(D a g) = lam(g) + up for all g in gs"""
    for _ in range(tries):
        a = letter(rng)
        if all(lam(mv(a, g)) == lam(g) + up for g in gs): return a
    return None
def to_K(rng, g, tries=4000):
    """letter mu with lag(D mu g) = LK (g at layer 3 adjacent to LK)"""
    for _ in range(tries):
        a = letter(rng)
        if lag(mv(a, g)) == LK: return a
    return None
def dist(L1, L2):
    return 4 - (len(L1 & L2).bit_length() - 1)

# ---- exact lifts (for the closure / Groebner stage)
def letter_full(rng):
    """random engine letter with its data: (S, P, Bm) ; P = list of 4 column images, Bm symmetric 4x4"""
    P = gl4(rng)
    Bm = [[0] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(i, 4):
            Bm[i][j] = Bm[j][i] = rng.randrange(2)
    return (symp_from(P, Bm), P, Bm)
def symp_from(P, Bm):
    Pm = {}
    for v in range(16):
        w = 0
        for j in range(4):
            if (v >> j) & 1: w ^= P[j]
        Pm[v] = w
    Pinv = {w: v for v, w in Pm.items()}
    PT = [sum(((Pinv[1 << j] >> i) & 1) << j for j in range(4)) for i in range(4)]
    def Bx(x):
        r = 0
        for i in range(4):
            if sum(Bm[i][j] * ((x >> j) & 1) for j in range(4)) & 1: r |= 1 << i
        return r
    S = [P[i] | (Bx(P[i]) << 4) for i in range(4)] + [PT[i] << 4 for i in range(4)]
    return tuple(S)
def to_matrix(P, Bm, shift=0, zph=0):
    """exact 16x16 monomial: e_x -> i^q(Px) e_(Px + shift) * (-1)^(zph . Px); q(y) = sum B_jj y_j + 2 sum_{j<k} B_jk y_j y_k"""
    import numpy as np
    M = np.zeros((16, 16), dtype=complex)
    for x in range(16):
        y = 0
        for j in range(4):
            if (x >> j) & 1: y ^= P[j]
        yb = [(y >> j) & 1 for j in range(4)]
        q = sum(Bm[j][j] * yb[j] for j in range(4)) + 2 * sum(Bm[j][k] * yb[j] * yb[k] for j in range(4) for k in range(j + 1, 4))
        q += 2 * (bin(zph & y).count('1') & 1)
        M[y ^ shift, x] = [1, 1j, -1, -1j][q % 4]
    return M
