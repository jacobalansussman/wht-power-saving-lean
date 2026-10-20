"""t_face3.py (rm-theory): an independent recomputation of the block-3 FACE certificate (wht38/tiny/close), in the form
needed for Theorem C3 (NOTE 8.9): f = 3 digit columns, restricted to the antisymmetric part Lambda^3(C^8) (kept by every
f-fold letter, since a letter acts the same way on every column), N = span of a(T) over the 24 face triples T.

Computes, exactly (ranks modulo two primes p = 1 mod 4 with i -> sqrt(-1), and over Q(i) with sympy for the 7 caps and
3 of the demands as a cross-check):
  cap(w)  = dim(N + dir(w) N) - dim N  for the 7 directions;
  dem(b)  = dim(N + K_b N) - dim N     with K_b = dir(b1) dir(b2) dir(b3), for all 28 bases b of F_2^3;
  leak(h) for the 168 address maps (sanity check against theory-guide/out_leak_m3.txt: {0: 6, 8: 54, 16: 90, 24: 18});
  that N is kept by the 7 shifts and by random diagonal matrices.
dir(w) = ((1+i) I + (1-i) X^w) / 2; we use 2 dir(w) (same ranks).
"""
import itertools, random
import numpy as np

n = 8
TRI = list(itertools.combinations(range(n), 3))
TIDX = {T: k for k, T in enumerate(TRI)}
FACES = [T for T in TRI if any(len({(t >> b) & 1 for t in T}) == 1 for b in range(3))]
PRIMES = [2147483497, 1000000009]          # both prime, = 1 mod 4


def sqrtm1(p):
    for a in range(2, 10000):
        r = pow(a, (p - 1) // 4, p)
        if r * r % p == p - 1: return r


def compound3(Mf, p):
    """Lambda^3 of an 8x8 matrix over F_p (entries ints mod p): 56x56 matrix of 3x3 minors, column T' -> row T."""
    C = np.zeros((56, 56), dtype=object)
    for a, T in enumerate(TRI):
        for b, U in enumerate(TRI):
            sub = [[Mf[t][u] for u in U] for t in T]
            d = (sub[0][0] * (sub[1][1] * sub[2][2] - sub[1][2] * sub[2][1])
                 - sub[0][1] * (sub[1][0] * sub[2][2] - sub[1][2] * sub[2][0])
                 + sub[0][2] * (sub[1][0] * sub[2][1] - sub[1][1] * sub[2][0])) % p
            C[a, b] = d
    return C


def rank_mod(A, p):
    A = [list(map(int, row)) for row in A]
    r = 0; rows = len(A); cols = len(A[0])
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] % p), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c]:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        r += 1
        if r == rows: break
    return r


def dirmat(w, p, I):
    M = [[0] * n for _ in range(n)]
    for x in range(n):
        M[x][x] = (M[x][x] + 1 + I) % p
        M[x ^ w][x] = (M[x ^ w][x] + 1 - I) % p
    return M


def matmul(A, B, p):
    return [[sum(A[i][k] * B[k][j] for k in range(n)) % p for j in range(n)] for i in range(n)]


def leak(C, p):
    """f_N(g) for g with compound C: rank of [N | C N] - dim N, N = face columns."""
    Ncols = [TIDX[T] for T in FACES]
    Nmat = np.zeros((56, len(Ncols)), dtype=object)
    for j, k in enumerate(Ncols): Nmat[k, j] = 1
    CN = C[:, Ncols] % p
    big = np.concatenate([Nmat, CN], axis=1)
    return rank_mod(big.T.tolist(), p) - len(Ncols)


def run():
    print(f"face triples: {len(FACES)} of {len(TRI)}")
    res = {}
    for p in PRIMES:
        I = sqrtm1(p)
        caps = {w: leak(compound3(dirmat(w, p, I), p), p) for w in range(1, n)}
        bases = []
        for cols in itertools.combinations(range(1, n), 3):
            if rank2(list(cols)) == 3: bases.append(cols)
        dems = {}
        for b in bases:
            K = dirmat(b[0], p, I)
            for w in b[1:]: K = matmul(dirmat(w, p, I), K, p)
            dems[b] = leak(compound3(K, p), p)
        # address maps
        leaks = {}
        for cols in itertools.product(range(1, n), repeat=3):
            if rank2(list(cols)) != 3: continue
            P = [[0] * n for _ in range(n)]
            for x in range(n):
                y = 0
                for k in range(3):
                    if (x >> k) & 1: y ^= cols[k]
                P[y][x] = 1
            lk = leak(compound3(P, p), p); leaks[lk] = leaks.get(lk, 0) + 1
        # shifts and random diagonals keep N
        keep = all(leak(compound3([[1 if y == x ^ a else 0 for x in range(n)] for y in range(n)], p), p) == 0 for a in range(1, n))
        rng = random.Random(1)
        diag_ok = True
        for _ in range(5):
            D = [[(rng.randrange(1, p) if x == y else 0) for x in range(n)] for y in range(n)]
            diag_ok &= leak(compound3(D, p), p) == 0
        res[p] = (caps, dems, leaks, keep, diag_ok)
        print(f"p = {p}: caps {caps}; bases {len(bases)}: dem min {min(dems.values())} max {max(dems.values())}; "
              f"address-map leaks {dict(sorted(leaks.items()))}; shifts keep N: {keep}; diagonals keep N: {diag_ok}",
              flush=True)
    a, b = res[PRIMES[0]], res[PRIMES[1]]
    print("two primes agree:", a[0] == b[0] and a[1] == b[1] and a[2] == b[2])


def rank2(v):
    basis = []
    for x in v:
        for bb in basis: x = min(x, x ^ bb)
        if x: basis.append(x)
    return len(basis)


if __name__ == "__main__":
    run()
