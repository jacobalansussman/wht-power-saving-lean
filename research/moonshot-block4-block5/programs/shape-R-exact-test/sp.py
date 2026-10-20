"""sp.py -- o3-close: exact symplectic model of the block-4 engine game (Clifford group mod Paulis -> Sp(8,2)).

Vectors of V = F_2^8 are ints: bits 0..3 = X-part x, bits 4..7 = Z-part z.  Form <(x,z),(x',z')> = x.z' + z.x'.
A term h (product of engine letters and D = dir(e1) = sqrt(X_e1)) acts on Paulis by conjugation; sigma(h) in Sp(8,2).
  L0 = Z-type Lagrangian (computational basis).  Free letters B' = monomial Cliffords = Stab(L0) -> parabolic P0.
  D = dir(e1) -> transvection tau by v = (e1|0).
  STATE of a term h: (mu, kappa) = (sigma(h) L0, sigma(h) Lambda_K), Lambda_K = sigma(K)^-1 L0 (always opposite).
  depth(h) = dist(mu, L0); distance to the kernel = dist(kappa, L0).  Left free letter b: state -> sigma(b) state.
  D: state -> tau state.  (Checked against the 16x16 matrices of bm.py in check_model.py.)
Lagrangians are frozensets of 16 ints.  dist(A, B) = 4 - dim(A cap B)."""
import itertools, random

M = 4
NV = 1 << (2 * M)

def form(a, b):
    x, z = a & 15, a >> 4
    x2, z2 = b & 15, b >> 4
    return (bin(x & z2).count('1') + bin(z & x2).count('1')) & 1

def span(vs):
    S = {0}
    for v in vs:
        if v not in S:
            S |= {s ^ v for s in S}
    return frozenset(S)

def dim(S):
    return len(S).bit_length() - 1

def dist(A, B):
    return M - dim(A & B)

L0 = span([1 << (4 + j) for j in range(M)])          # Z-type

def transvection(v):
    return lambda w: w ^ v if form(w, v) else w

def apply_map(f, S):
    return frozenset(f(w) for w in S)

# linear maps as tuples of images of the 8 basis vectors
def mat_apply(Mt, w):
    out = 0
    j = 0
    while w:
        if w & 1:
            out ^= Mt[j]
        w >>= 1; j += 1
    return out

def mat_mul(A, B):        # (A B)(w) = A(B(w))
    return tuple(mat_apply(A, B[j]) for j in range(2 * M))

def mat_inv(A):
    # brute force via solving: build inverse by Gaussian elimination over F2 on 8x8
    n = 2 * M
    rows = [(A[j], 1 << j) for j in range(n)]   # image, preimage-basis
    # we want B with B(A e_j) = e_j; columns: solve
    basis = []
    for img, pre in rows:
        for (bi, bp, piv) in basis:
            if img >> piv & 1:
                img ^= bi; pre ^= bp
        assert img != 0
        piv = img.bit_length() - 1
        basis = [(bi ^ img, bp ^ pre, p) if bi >> piv & 1 else (bi, bp, p) for (bi, bp, p) in basis]
        basis.append((img, pre, piv))
    inv = [0] * n
    for (bi, bp, piv) in basis:
        assert bi == 1 << piv
        inv[piv] = bp
    return tuple(inv)

def mat_of(f):
    return tuple(f(1 << j) for j in range(2 * M))

def is_symplectic(A):
    for i in range(2 * M):
        for j in range(2 * M):
            if form(A[i], A[j]) != form(1 << i, 1 << j):
                return False
    return True

TAU = mat_of(transvection(1))                 # D = dir(e1): v = (e1|0)
def tau_z(z):
    return mat_of(transvection(z))            # dir(z), z in X-part bits 0..3

# address map A in GL(4,2) (columns A[j] = image of e_j, ints 0..15): x -> A x, z -> A^-T z
def gl_cols_to_sp(Acols):
    # A^-T: compute inverse of A as 4x4 then transpose
    n = M
    Ainv = [0] * n
    # inverse via brute force over 4x4
    def apply4(cols, x):
        o = 0
        for j in range(n):
            if x >> j & 1:
                o ^= cols[j]
        return o
    for y in range(16):
        pass
    inv_map = {}
    for x in range(16):
        inv_map[apply4(Acols, x)] = x
    Ainv = [inv_map[1 << j] for j in range(n)]
    # A^-T e_j : (A^-T)_{ij} = (A^-1)_{ji}; column j of A^-T = row j of A^-1
    AinvT = []
    for j in range(n):
        col = 0
        for i in range(n):
            # (A^-1)_{j i} = bit j of Ainv[i]
            if Ainv[i] >> j & 1:
                col |= 1 << i
        AinvT.append(col)
    imgs = [Acols[j] for j in range(n)] + [AinvT[j] << 4 for j in range(n)]
    return tuple(imgs)

def sym_to_sp(S):
    """phase letter with symmetric matrix S (list of 4 row ints): (x, z) -> (x, z + S x)"""
    imgs = []
    for j in range(M):
        imgs.append((1 << j) | (S[j] << 4))     # column j of S (S symmetric: row j = column j)
    for j in range(M):
        imgs.append(1 << (4 + j))
    return tuple(imgs)

def rand_gl(rng):
    while True:
        cols = [rng.randrange(1, 16) for _ in range(M)]
        if len(span(cols)) == 16:
            return cols

def rand_sym(rng):
    S = [[0] * M for _ in range(M)]
    for i in range(M):
        for j in range(i, M):
            b = rng.randrange(2)
            S[i][j] = S[j][i] = b
    return [sum(S[i][j] << j for j in range(M)) for i in range(M)]

def rand_P0(rng):
    return mat_mul(gl_cols_to_sp(rand_gl(rng)), sym_to_sp(rand_sym(rng)))

ID = tuple(1 << j for j in range(2 * M))
KSP = ID
for z in (1, 2, 4, 8):
    KSP = mat_mul(tau_z(z), KSP)
LK = apply_map(lambda w: mat_apply(mat_inv(KSP), w), L0)     # Lambda_K = sigma(K)^-1 L0
LD = apply_map(lambda w: mat_apply(TAU, w), L0)              # Lambda_D = tau L0

def img(A, S):
    return frozenset(mat_apply(A, w) for w in S)

# ---- all Lagrangians (orbit of L0 under Sp(8,2) generated by P0 generators and tau)
def all_lagrangians(seed=1):
    rng = random.Random(seed)
    gens = [TAU] + [rand_P0(rng) for _ in range(6)]
    seen = {L0}; frontier = [L0]
    while frontier:
        nf = []
        for L in frontier:
            for g in gens:
                L2 = img(g, L)
                if L2 not in seen:
                    seen.add(L2); nf.append(L2)
        frontier = nf
    return seen

def neighbours_of_L0(LAG):
    return [L for L in LAG if dist(L, L0) == 1]

if __name__ == '__main__':
    assert is_symplectic(TAU)
    rng = random.Random(3)
    for _ in range(20):
        P = rand_P0(rng); assert is_symplectic(P); assert img(P, L0) == L0
    LAG = all_lagrangians()
    print("Lagrangians:", len(LAG))
    from collections import Counter
    print("distance distribution from L0:", sorted(Counter(dist(L, L0) for L in LAG).items()))
    print("dist(L0, Lambda_K) =", dist(L0, LK), " dist(L0, Lambda_D) =", dist(L0, LD))
    nb = neighbours_of_L0(LAG)
    print("neighbours of L0:", len(nb))
