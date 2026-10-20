"""sp2.py (w3-b5): exact F_2 symplectic toolkit for block m (default 5): Sp(2m, 2), Lagrangians (frames), the dual polar
graph, the images of the engine letters.  Pure Python (no numpy), so it runs on the Mac.

Convention (matches mw-big56/mwlib.py): a Pauli X(a)Z(b) <-> the integer v = a | (b << m) (bits 0..m-1 = X part on
qubits 0..m-1, bits m..2m-1 = Z part).  omega(u, v) = a_u.b_v + b_u.a_v (mod 2).  A Clifford g acts by conjugation
P -> g P g^-1; gbar is that map on F_2^{2m}.  L0 = Z-type Lagrangian = frame of the identity; the frame of a group element g
is F(g) = gbar^-1 L0 (left B-cosets <-> frames), and layer(g) = d(L0, F(g)) = m - dim(L0 cap F(g)).
D = dir(e1) = sqrt(X) on qubit 0: Z_0 -> X_0 Z_0 (a transvection).  K = product of all dir(e_j): Z_j -> X_j Z_j.
The free letters B map onto P = Stab(L0) = {(a, b) -> (A a, A^-T b + S A a)}, A in GL(m,2), S symmetric."""
import random

M = 5
NB = 2 * M
MASKX = (1 << M) - 1


def par(x):
    return bin(x).count('1') & 1


def omega(u, v):
    return par((u & MASKX) & (v >> M)) ^ par((u >> M) & (v & MASKX))


# ---------------------------------------------------------------- linear maps on F_2^{2m}: list of images of basis vectors
def app(Mp, v):
    y = 0; i = 0
    while v:
        if v & 1:
            y ^= Mp[i]
        v >>= 1; i += 1
    return y


def mul(A, B):           # (A B)(v) = A(B(v))
    return [app(A, b) for b in B]


IDM = [1 << i for i in range(NB)]


def inv(A):
    # Gaussian elimination on the augmented system: find X with A X = I.  Columns of A are A[i].
    n = NB
    rows = []      # represent the matrix by rows: row r = sum over columns i of bit r of A[i]
    for r in range(n):
        row = 0
        for i in range(n):
            if (A[i] >> r) & 1:
                row |= 1 << i
        rows.append([row, 1 << r])
    piv = 0
    for c in range(n):
        pr = next((r for r in range(piv, n) if (rows[r][0] >> c) & 1), None)
        assert pr is not None, "singular"
        rows[piv], rows[pr] = rows[pr], rows[piv]
        for r in range(n):
            if r != piv and (rows[r][0] >> c) & 1:
                rows[r][0] ^= rows[piv][0]; rows[r][1] ^= rows[piv][1]
        piv += 1
    # now rows[c] = [e_c, (A^-1 row c)]; A^-1 column i = bits i of rows
    out = []
    for i in range(n):
        col = 0
        for r in range(n):
            if (rows[r][1] >> i) & 1:
                col |= 1 << r
        out.append(col)
    return out


def is_symp(A):
    for i in range(NB):
        for j in range(NB):
            if omega(A[i], A[j]) != omega(1 << i, 1 << j):
                return False
    return True


def order(A, cap=100000):
    R = list(A); k = 1
    while R != IDM:
        R = mul(A, R); k += 1
        if k > cap:
            return None
    return k


# ---------------------------------------------------------------- subspaces / Lagrangians
def rref(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            if v ^ b < v:
                v ^= b
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    # full reduction
    basis.sort(reverse=True)
    for i in range(len(basis)):
        hb = basis[i].bit_length() - 1
        for j in range(len(basis)):
            if j != i and (basis[j] >> hb) & 1:
                basis[j] ^= basis[i]
    return tuple(sorted(basis, reverse=True))


def rank(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            if v ^ b < v:
                v ^= b
        if v:
            basis.append(v); basis.sort(reverse=True)
    return len(basis)


def dist(L1, L2):
    return rank(list(L1) + list(L2)) - M


def image(A, L):
    return rref([app(A, v) for v in L])


def span_all(L):
    out = {0}
    for v in L:
        out |= {x ^ v for x in out}
    return out


L0 = rref([1 << (M + j) for j in range(M)])
D = list(IDM); D[M] = (1 << M) | 1                       # Z_0 -> X_0 Z_0
K = list(IDM)
for j in range(M):
    K[M + j] = (1 << (M + j)) | (1 << j)
Fstar = image(inv(K), L0)                                # F(K) = Kbar^-1 L0


def neighbors(L):
    Ls = span_all(L); out = set()
    for w in range(1, 1 << NB):
        if w in Ls:
            continue
        H = [x for x in L if omega(x, w) == 0]
        # L cap w^perp: solve via span (L has 32 elements)
        Hs = [x for x in Ls if omega(x, w) == 0]
        out.add(rref(Hs + [w]))
    return out


# ---------------------------------------------------------------- the parabolic P = Stab(L0)
def levi(A5):
    """A5: list of m images (m-bit ints) of the X basis vectors (an invertible m x m matrix).  (a, b) -> (A a, A^-T b)."""
    # A^-T: compute the inverse of A5 as an m x m matrix, then transpose
    n = M
    rows = []
    for r in range(n):
        row = 0
        for i in range(n):
            if (A5[i] >> r) & 1:
                row |= 1 << i
        rows.append([row, 1 << r])
    piv = 0
    for c in range(n):
        pr = next((r for r in range(piv, n) if (rows[r][0] >> c) & 1), None)
        if pr is None:
            return None
        rows[piv], rows[pr] = rows[pr], rows[piv]
        for r in range(n):
            if r != piv and (rows[r][0] >> c) & 1:
                rows[r][0] ^= rows[piv][0]; rows[r][1] ^= rows[piv][1]
        piv += 1
    Ainv_rows = [rows[r][1] for r in range(n)]          # row r of A^-1 as a bitmask over columns
    # A^-T column i = row i of A^-1
    out = [A5[i] for i in range(n)] + [Ainv_rows[i] << M for i in range(n)]
    return out


def unip(S):
    """S: symmetric m x m over F2 as list of row bitmasks.  (a, b) -> (a, b + S a)."""
    out = []
    for i in range(M):
        out.append((1 << i) | (S[i] << M))
    for i in range(M):
        out.append(1 << (M + i))
    return out


def rand_gl(rng):
    while True:
        A = [rng.randrange(1, 1 << M) for _ in range(M)]
        if rank(A) == M:
            return A


def rand_sym(rng):
    S = [0] * M
    for i in range(M):
        for j in range(i, M):
            if rng.random() < 0.5:
                S[i] |= 1 << j; S[j] |= 1 << i
    return S


def rand_P(rng):
    return mul(levi(rand_gl(rng)), unip(rand_sym(rng)))


def selftest():
    rng = random.Random(1)
    assert is_symp(D) and is_symp(K)
    assert dist(L0, image(D, L0)) == 1 and dist(L0, Fstar) == M
    for _ in range(20):
        Pm = rand_P(rng)
        assert is_symp(Pm) and image(Pm, L0) == L0
        assert mul(Pm, inv(Pm)) == IDM
    nb = neighbors(L0)
    assert len(nb) == 2 * (2 ** M - 1), len(nb)
    return True


if __name__ == "__main__":
    print("selftest", selftest(), "neighbors of L0:", len(neighbors(L0)))
