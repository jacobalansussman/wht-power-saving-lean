# b56-ceilings: the K-invariant piece of the hook certificate (EXACT check over Z[i], small cases).
# Hook: W = V (x) Lambda^k V, N = span{e_x (x) e_T : x in T}.  Claim [proved by the identity below]:
#   psi(w) = sum_x (-1)^{|x|} e_x (x) (e_x ^ w),  w in Lambda^(k-1),  lies in N, and
#   (K (x) K^k) psi(w) = i^m psi(K^(k-1) w)      (because K Z_1 K^T = i^m Z_1, K symmetric, Z_1 = diag((-1)^|x|)),
# so psi(Lambda^(k-1)) (dimension C(n,k-1)) is a K-invariant subspace of N, N cap KN contains it, and
#   dem(hook, k) <= k C(n,k) - C(n,k-1).
# This file checks the identity K Z_1 K = i^m Z_1 exactly (Gaussian integers, 2^m K), the injectivity of psi, and the
# K-invariance of psi(Lambda^(k-1)) inside V (x) Lambda^k exactly for the given (m, k).
# Usage: python3 -I hook_kinv.py M K
import sys, itertools
m, k = int(sys.argv[1]), int(sys.argv[2]); n = 2 ** m
def gm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gadd(a, b): return (a[0] + b[0], a[1] + b[1])
# 2 dir(z) = (1+i) I + (1-i) X_z ; 2^m K = prod_j 2 dir(e_j): entry (y, x) = prod_j ((1+i) if bit j equal else (1-i))
def K2(y, x):
    c = (1, 0)
    for j in range(m):
        c = gm(c, (1, 1) if ((x ^ y) >> j) & 1 == 0 else (1, -1))
    return c
Kmat = [[K2(y, x) for x in range(n)] for y in range(n)]
sym = all(Kmat[y][x] == Kmat[x][y] for x in range(n) for y in range(n))
z = [(-1) ** bin(x).count('1') for x in range(n)]
# (2^m K) Z_1 (2^m K) = 4^m i^m Z_1 ?
im = [(1, 0), (0, 1), (-1, 0), (0, -1)][m % 4]
ok = True
for y in range(n):
    for y2 in range(n):
        s = (0, 0)
        for x in range(n):
            s = gadd(s, gm(gm(Kmat[y][x], (z[x], 0)), Kmat[x][y2]))
        want = gm((4 ** m * (z[y] if y == y2 else 0), 0), im)
        ok &= (s == want)
print(f"m={m}: K symmetric {sym};  K Z_1 K = i^m Z_1 exactly: {ok}")
# psi(w) for basis w = e_U, U a (k-1)-subset: sparse dict (x, T) -> sign * z_x
subs_km1 = list(itertools.combinations(range(n), k - 1))
def wedge_sign(x, U):
    # e_x ^ e_U = sign * e_{sorted(U+x)}
    return -1 if sum(1 for u in U if u < x) % 2 else 1
def psi(U):
    v = {}
    for x in range(n):
        if x in U: continue
        T = tuple(sorted(U + (x,)))
        v[(x, T)] = z[x] * wedge_sign(x, U)
    return v
inN = all(x in T for U in subs_km1 for (x, T) in psi(U))
# injectivity: distinct U give disjoint supports? (x, T) determines U = T \ {x}
supports = [frozenset(psi(U)) for U in subs_km1]
inj = all(len(s) > 0 for s in supports) and len(set().union(*supports)) == sum(len(s) for s in supports)
print(f"m={m} k={k}: psi(Lambda^(k-1)) inside N: {inN};  psi injective (disjoint nonempty supports): {inj};  dim = C(n,k-1) = {len(subs_km1)}")
