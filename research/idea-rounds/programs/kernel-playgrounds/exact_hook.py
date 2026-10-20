# playgrounds/exact_hook.py -- EXACT (over Q) cap and demand of the f=5 hook at block 3 for integer kernels.
# Same construction as hookk.py, but the wedge matrices and ranks are computed with exact integers
# (Bareiss fraction-free elimination), so cap and dem are the true ranks over Q (= over C for rational matrices).
# Usage: python3 -I exact_hook.py
import itertools, sys
from fractions import Fraction

def kron(A, B):
    return [[a * b for a in ra for b in rb] for ra in A for rb in B]

def det_int(M):
    # Bareiss determinant for small integer matrices
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None: return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def wedge(g, k, subsets):
    return [[det_int([[g[t][u] for u in U] for t in T]) for U in subsets] for T in subsets]

def rank_int(A):
    # Bareiss-style fraction-free row reduction, exact rank over Q
    A = [row[:] for row in A]
    rows = len(A); cols = len(A[0]) if rows else 0
    r = 0; prev = 1
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        pr = A[r]; pv = pr[c]
        for i in range(r + 1, rows):
            ai = A[i]; aic = ai[c]
            if aic == 0:
                A[i] = [(x * pv) // prev for x in ai] if prev != 1 else [x * pv for x in ai]
                continue
            A[i] = [(ai[j] * pv - aic * pr[j]) // prev for j in range(cols)]
        prev = pv; r += 1
        if r == rows: break
    return r

def hook_exact(K, m=3, k=4):
    n = 2 ** m
    I2 = [[1, 0], [0, 1]]
    K1 = K
    for _ in range(m - 1): K1 = kron(K1, I2)
    T = K
    for _ in range(m - 1): T = kron(T, K)
    subsets = list(itertools.combinations(range(n), k))
    basis = [(x, S) for x in range(n) for S in subsets]
    inN = [x in S for (x, S) in basis]
    Nidx = [i for i, b in enumerate(inN) if b]; Oidx = [i for i, b in enumerate(inN) if not b]
    out = {}
    for name, g in (('cap', K1), ('dem', T)):
        Wg = kron(g, wedge(g, k, subsets))
        B = [[Wg[i][j] for j in Nidx] for i in Oidx]
        out[name] = rank_int(B)
    return out

if __name__ == '__main__':
    for name, K in (('H', [[1, 1], [1, -1]]), ('L', [[1, 0], [1, 1]]), ('Kinv', [[1, 0], [1, -1]])):
        r = hook_exact(K)
        print(name, 'exact over Q: cap =', r['cap'], ' dem =', r['dem'], ' ratio =', Fraction(r['dem'], r['cap']), flush=True)
