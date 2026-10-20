# indep_venn_check.py (meta/literature, 14 Oct 2026).  Independent check (own code, no import from venn_mask.py) of
# f-values of a two-factor Venn certificate at block m = 3:  X = Lambda^a(V) (x) Lambda^b(V), V = Q^8 (address bits
# 0,1,2), N = span{ e_T1 (x) e_T2 : |T1 cap T2| in J } (a coordinate subspace).  For a matrix g on V,
#   f(g) = dim(N + gN) - dim N = rank of the block of Lambda^a(g) (x) Lambda^b(g) with rows OUTSIDE N, columns IN N.
# Letters: H_j = the unnormalised 2x2 Hadamard on address bit j (scalars never change f), so
#   cap = f(H_0), f(rank-2 target, one spare row) = f(H_0 H_1), f(kernel) = f(H_0 H_1 H_2).
# Ranks mod two primes near 2^31 (a mod-p rank is <= the rational rank; agreement of two primes + the closed-form exact
# cap is the check).  Usage: python3 -I indep_venn_check.py a b J   (J like 0,2)
import sys, itertools
import numpy as np
a, b = int(sys.argv[1]), int(sys.argv[2]); J = {int(x) for x in sys.argv[3].split(',')}
n = 8
def hadamard_on(bits):
    g = np.eye(n, dtype=np.int64)
    for j in bits:
        h = np.zeros((n, n), dtype=np.int64)
        for x in range(n):
            y = x ^ (1 << j)
            h[x, x] = -1 if (x >> j) & 1 else 1
            h[y, x] = 1
        g = h @ g
    return g
def compound(g, k):
    S = list(itertools.combinations(range(n), k))
    C = np.zeros((len(S), len(S)), dtype=np.int64)
    for i, R in enumerate(S):
        sub = g[list(R)]
        for j, Cc in enumerate(S):
            C[i, j] = int(round(np.linalg.det(sub[:, list(Cc)].astype(float))))
    return S, C
def rank_mod(A, p):
    A = (A % p).astype(np.int64); r = 0; rows, cols = A.shape
    for c in range(cols):
        if r == rows: break
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0: continue
        pr = r + piv[0]
        if pr != r: A[[r, pr]] = A[[pr, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        nz = np.nonzero(A[:, c])[0]; nz = nz[nz != r]
        if len(nz): A[nz] = (A[nz] - (A[nz, c][:, None] * A[r][None, :]) % p) % p
        r += 1
    return r
Sa = list(itertools.combinations(range(n), a)); Sb = list(itertools.combinations(range(n), b))
pairs = [(i, j) for i in range(len(Sa)) for j in range(len(Sb))]
inN = np.array([len(set(Sa[i]) & set(Sb[j])) in J for i, j in pairs])
print(f"a={a} b={b} J={sorted(J)}: dim X {len(pairs)}, dim N {int(inN.sum())}", flush=True)
for name, bits in (("cap = f(H0)", [0]), ("f(H0 H1) rank-2 target", [0, 1]), ("f(H0 H1 H2) kernel", [0, 1, 2])):
    g = hadamard_on(bits)
    _, Ca = compound(g, a); _, Cb = compound(g, b)
    big = np.kron(Ca, Cb)
    blk = big[np.ix_(~inN, inN)]
    rks = [rank_mod(blk, p) for p in (2147483629, 2147483587)]
    print(f"  {name}: ranks mod two primes {rks}", flush=True)
