# v-ut-direct-logic e2_counter.py (pure python, F2 matrices of size 2mf <= 96; tiny, < 10 MB)
# Tests PROOF.md 7.2 (E2) "Q is a group" for the definition used in 7.1 and in symbol_check.py check (a):
#   M(k_f) - T_f(M_k) vanishes on deep rows x deep columns (deviation supported in boundary block rows OR columns).
# Letters (m = 2, address maps on row 0 only; label index = 2m*j + (0..m-1 for x, m..2m-1 for z), as symbol_check.py):
#   k : x_0[0] ^= x_1[0] + ... + x_{f-1}[0]   (address map A = I + e0 u^T)
#   l : x_j[0] ^= x_0[0] for all j >= 1       (address map B = I + u e0^T)
# Both are monomial Clifford address maps (in AGL(mf, F2)); both are in Q with r = 1 and symbol I.
# Claim tested: k*l is NOT in Q for any r (deep x deep deviation), but rank(M(kl) - I) stays bounded (repair).
m = 2
def ident(n): return [[int(i == j) for j in range(n)] for i in range(n)]
def mul(A, B):
    n, p = len(A), len(B[0]); Bt = list(zip(*B))
    return [[sum(a & b for a, b in zip(A[i], Bt[j])) & 1 for j in range(p)] for i in range(n)]
def T(A): return [list(r) for r in zip(*A)]
def inv(M):
    n = len(M); A = [M[i][:] + ident(n)[i] for i in range(n)]; r = 0
    for c in range(n):
        p = next(i for i in range(r, n) if A[i][c]); A[r], A[p] = A[p], A[r]
        for i in range(n):
            if i != r and A[i][c]: A[i] = [x ^ y for x, y in zip(A[i], A[r])]
        r += 1
    return [row[n:] for row in A]
def rank(M):
    A = [r[:] for r in M]; rk = 0; n = len(A); cols = len(A[0])
    for c in range(cols):
        p = next((i for i in range(rk, n) if A[i][c]), None)
        if p is None: continue
        A[rk], A[p] = A[p], A[rk]
        for i in range(n):
            if i != rk and A[i][c]: A[i] = [x ^ y for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk
def label(X):  # address map X (mf x mf) -> label matrix on (x|z): X on x-labels, X^-T on z-labels
    f = len(X) // m; Zt = T(inv(X)); M = [[0] * (2 * m * f) for _ in range(2 * m * f)]
    for j in range(f):
        for k in range(f):
            for a in range(m):
                for b in range(m):
                    M[2*m*j + a][2*m*k + b] = X[m*j + a][m*k + b]
                    M[2*m*j + m + a][2*m*k + m + b] = Zt[m*j + a][m*k + b]
    return M
def symplectic(M, f):
    J = [[0] * (2*m*f) for _ in range(2*m*f)]
    for j in range(f):
        for i in range(m): J[2*m*j + i][2*m*j + m + i] = 1; J[2*m*j + m + i][2*m*j + i] = 1
    return mul(mul(T(M), J), M) == J
def deep_dev(M, f, r):  # nonzero entries of M - I with BOTH block row and block column in [r, f - r)
    n = 2 * m * f
    return [(i, j) for i in range(n) for j in range(n) if (M[i][j] ^ (i == j)) and r <= i // (2*m) < f - r and r <= j // (2*m) < f - r]
ok = True
for f in (8, 12):
    n = m * f
    A = ident(n); B = ident(n)
    for j in range(1, f): A[0][m * j] = 1; B[m * j][0] = 1   # row-0 bit of column j is index m*j
    Mk, Ml = label(A), label(B); Mkl = mul(Mk, Ml)
    s = symplectic(Mk, f) and symplectic(Ml, f)
    k_in = not deep_dev(Mk, f, 1); l_in = not deep_dev(Ml, f, 1)
    kl_dev = {r: len(deep_dev(Mkl, f, r)) for r in range(1, f // 2)}
    rk = rank([[Mkl[i][j] ^ (i == j) for j in range(2*n)] for i in range(2*n)])
    print(f"f={f}: k,l symplectic: {s}; k in Q(r=1, symbol I): {k_in}; l in Q(r=1, symbol I): {l_in}")
    print(f"   k*l: deep x deep nonzero entries of M - I for r = 1..{f//2 - 1}: {kl_dev}  (all > 0 => k*l not in Q for any r)")
    print(f"   k*l: rank(M(kl) - I) = {rk}  (bounded: the bounded-rank repair of Q keeps symbol I)")
    ok &= s and k_in and l_in and all(v > 0 for v in kl_dev.values()) and rk <= 4
print("COUNTEREXAMPLE TO E2 AS WRITTEN CONFIRMED" if ok else "NOT CONFIRMED")
