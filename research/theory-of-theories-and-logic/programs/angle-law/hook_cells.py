# Skeptic's own hook computation (independent code).
# (A) Theorem H check: phi(g) = ||(1-P)GP||_F^2 vs C(n-2,k-1)(n - sum|g|^4) for random unitaries, small (n,k).
# (B) Cell reduction: for g = h^{(x)j} (x) I on n = 2^m points, the depth-j hook rank and phi at block m are
#     2^(m-j) * sum_a r_j(a) * C(n - 2^j, k - a), where r_j(a) is the FULL-depth hook rank at block j with k = a.
#     r_j(a) computed by float SVD (with the singular-value gap printed) and mod two primes p = 1 mod 4.
import itertools, sys
import numpy as np
from math import comb

def wedge(g, k):
    n = g.shape[0]
    subs = list(itertools.combinations(range(n), k))
    W = np.empty((len(subs), len(subs)), dtype=g.dtype)
    for r, T in enumerate(subs):
        gT = g[list(T), :]
        for c, S in enumerate(subs):
            W[r, c] = np.linalg.det(gT[:, list(S)]) if k > 0 else 1
    return W, subs

def hook_block(g, k):
    n = g.shape[0]
    W, subs = wedge(g, k)
    G = np.kron(g, W)
    inN = np.array([x in T for x in range(n) for T in subs])
    return G[np.ix_(~inN, inN)], inN.sum()

def hdir(j):
    h = np.array([[(1+1j)/2, (1-1j)/2], [(1-1j)/2, (1+1j)/2]])
    g = np.array([[1.0+0j]])
    for _ in range(j): g = np.kron(g, h)
    return g

# ---------- (A) ----------
rng = np.random.default_rng(12345)
print("(A) Theorem H, random unitaries")
worst = 0
for (n, k) in [(3,1),(4,1),(4,2),(5,2),(6,2),(6,3),(7,3),(8,2),(8,4)]:
    for trial in range(3):
        Z = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
        Q, R = np.linalg.qr(Z)
        B, dimN = hook_block(Q, k)
        phi = np.linalg.norm(B)**2
        pred = comb(n-2, k-1) * (n - (np.abs(Q)**4).sum())
        worst = max(worst, abs(phi - pred))
    print(f"  n={n} k={k}: last phi={phi:.10f} pred={pred:.10f}")
print("  max |phi - closed form| over all trials:", worst)

# ---------- mod-p rank ----------
def modp_rank(M, p):
    M = M.copy() % p
    rows, cols = M.shape; r = 0
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        pr = r + nz[0]
        if pr != r: M[[r, pr]] = M[[pr, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p-2, p)) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): M[nzr] = (M[nzr] - (col[nzr, None] * M[r][None, :]) % p) % p
        r += 1
    return r

def det_modp(A, p):
    A = A.copy() % p; n = A.shape[0]; d = 1
    for c in range(n):
        nz = np.nonzero(A[c:, c])[0]
        if len(nz) == 0: return 0
        pr = c + nz[0]
        if pr != c: A[[c, pr]] = A[[pr, c]]; d = -d
        d = d * int(A[c, c]) % p
        inv = pow(int(A[c, c]), p-2, p)
        for r2 in range(c+1, n):
            if A[r2, c]:
                A[r2] = (A[r2] - A[r2, c] * inv % p * A[c]) % p
    return d % p

def hook_block_modp(j, k, p):
    I = next(pow(t, (p-1)//4, p) for t in range(2, 1000) if pow(t, (p-1)//2, p) == p-1)
    inv2 = (p+1)//2
    a_ = (1+I)*inv2 % p; b_ = (1-I)*inv2 % p
    h = np.array([[a_, b_], [b_, a_]], dtype=np.int64)
    g = np.array([[1]], dtype=np.int64)
    for _ in range(j): g = np.kron(g, h) % p
    n = g.shape[0]
    subs = list(itertools.combinations(range(n), k))
    W = np.zeros((len(subs), len(subs)), dtype=np.int64)
    for r, T in enumerate(subs):
        for c, S in enumerate(subs):
            W[r, c] = det_modp(g[np.ix_(list(T), list(S))], p) if k > 0 else 1
    G = np.kron(g, W) % p
    inN = np.array([x in T for x in range(n) for T in subs])
    return G[np.ix_(~inN, inN)]

# ---------- (B) ----------
print("\n(B) full-depth hook ranks r_j(a) at block j (float SVD gap, and mod two primes)")
PR = [1000033 if 1000033 % 4 == 1 else None, 999961]
PR = [p for p in [1000037, 999961, 1000033, 998989] if p % 4 == 1 and all(p % d for d in range(2, int(p**0.5)+1))][:2]
print("  primes:", PR)
r = {}; ph = {}
for j in (1, 2, 3):
    n = 2**j
    for a in range(0, n+1):
        B, dimN = hook_block(hdir(j), a)
        if B.size == 0:
            r[(j,a)] = 0; ph[(j,a)] = 0.0
            print(f"  j={j} a={a}: empty block (dimN={dimN})"); continue
        s = np.linalg.svd(B, compute_uv=False)
        rk = int((s > 1e-8).sum())
        gap_lo = s[rk-1] if rk > 0 else float('nan'); gap_hi = s[rk] if rk < len(s) else 0.0
        rks = []
        if B.shape[0]*B.shape[1] <= 400*400:
            for p in PR:
                rks.append(modp_rank(hook_block_modp(j, a, p), p))
        r[(j,a)] = rk; ph[(j,a)] = (s**2).sum()
        predphi = comb(n-2, a-1) * n * (1 - 2.0**-j) if a >= 1 else 0
        print(f"  j={j} a={a}: block {B.shape}, float rank {rk} (smallest kept sv {gap_lo:.3e}, largest dropped {gap_hi:.3e}), mod-p ranks {rks}; phi {ph[(j,a)]:.6f} (Thm H {predphi:.6f})")

def profile(m, k):
    n = 2**m
    out = []
    for j in range(1, min(m, 3)+1):
        rk = 2**(m-j) * sum(r[(j,a)] * comb(n - 2**j, k - a) for a in range(0, 2**j+1) if 0 <= k-a <= n-2**j)
        phi = 2**(m-j) * sum(ph[(j,a)] * comb(n - 2**j, k - a) for a in range(0, 2**j+1) if 0 <= k-a <= n-2**j)
        out.append((j, rk, phi))
    return out

print("\n(C) assembled profiles: (depth j, rank, phi) ")
for (m, k) in [(3,2),(3,4),(4,2),(4,5),(4,6),(4,7),(4,8),(4,10)]:
    pr = profile(m, k)
    s = ", ".join(f"d{j}: rank {rk} phi {phi:.3f} s={phi/rk:.6f}" for j, rk, phi in pr)
    print(f"  m={m} k={k}: {s}")
pr = profile(4, 8)
r1, r3 = pr[0][1], pr[2][1]
print(f"\nblock 4, k = 8: depth-3/depth-1 = {r3}/{r1} = {r3/r1:.6f}; needed > 2.5 i.e. rank_3 > {2.5*r1}; slack per array {2.5*r1 - r3}, two arrays {2*(2.5*r1 - r3)}")
print(f"5 moves x {r1} = {5*r1}; 2 x {r3} = {2*r3}; slack {5*r1 - 2*r3}")
