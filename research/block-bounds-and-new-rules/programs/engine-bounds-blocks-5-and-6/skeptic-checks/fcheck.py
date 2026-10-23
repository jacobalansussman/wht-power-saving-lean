# v-pf-induction-logic, 13 Oct 2026.  Independent check of the block-4 hook (k = 8) values
# f_1 = 27456, f_2 >= 48576, f_3 >= 68320 used as the base of PROOF.md Theorem 4.
# Method (independent of b56-compute's code): a rank-j block g = dir(e_1)...dir(e_j) acts as G (x) I on
# V = C^(2^j) (x) C^(2^(m-j)), so V = sum_c W_c (c = the untouched bits, W_c = C^d, d = 2^j) and g = sum_c G.
# Then Lambda^k V = sum over compositions (k_c) of tensor_c Lambda^(k_c) W_c, and the hook N = {e_x (x) e_T : x in T}
# splits as sum over (c0, k) of N1(W_c0, k_c0) (x) (full other factors).  Hence
#   f(g) = sum_{c0} sum_{k} phi_d(k_c0) * prod_{c != c0} C(d, k_c),
#   phi_d(j) = rank of the block of G (x) Lambda^j G with rows (y,S), y not in S, columns (x,T), x in T.
# G = a^j * B with B = [[1,-i],[-i,1]]^(x) j (a = (1+i)/2 is a scalar, rank unaffected); entries in Z[i].
# Ranks mod primes p = 1 mod 4 (i -> sqrt(-1) mod p): certified LOWER bounds on the char-0 rank; two primes.
import sys, itertools
from math import comb

def sqrtm1(p):
    for g in range(2, 1000):
        r = pow(g, (p - 1) // 4, p)
        if r * r % p == p - 1: return r
    raise ValueError

def rank_mod(M, p):
    M = [row[:] for row in M]; r = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(M)):
            if M[i][c] % p: piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [v * inv % p for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                t = M[i][c]; Mr = M[r]
                M[i] = [(vi - t * vr) % p for vi, vr in zip(M[i], Mr)]
        r += 1
    return r

def det_mod(A, p):
    A = [row[:] for row in A]; n = len(A); d = 1
    for c in range(n):
        piv = None
        for i in range(c, n):
            if A[i][c] % p: piv = i; break
        if piv is None: return 0
        if piv != c: A[c], A[piv] = A[piv], A[c]; d = -d
        d = d * A[c][c] % p; inv = pow(A[c][c], p - 2, p)
        for i in range(c + 1, n):
            if A[i][c] % p:
                t = A[i][c] * inv % p
                A[i] = [(vi - t * vc) % p for vi, vc in zip(A[i], A[c])]
    return d % p

def phi(j, kk, p, I):
    d = 2 ** j
    B1 = [[1, (-I) % p], [(-I) % p, 1]]
    G = [[1]]
    for _ in range(j):
        G = [[G[a][b] * B1[c][e] % p for b in range(len(G)) for e in range(2)] for a in range(len(G)) for c in range(2)]
    subs = list(itertools.combinations(range(d), kk))
    minors = {}
    for S in subs:
        for T in subs:
            minors[(S, T)] = det_mod([[G[s][t] for t in T] for s in S], p) if kk else 1
    rows = [(y, S) for S in subs for y in range(d) if y not in S]
    cols = [(x, T) for T in subs for x in T]
    if not rows or not cols: return 0
    M = [[G[y][x] * minors[(S, T)] % p for (x, T) in cols] for (y, S) in rows]
    return rank_mod(M, p)

def f_block(m, k, j, p):
    I = sqrtm1(p); d = 2 ** j; nb = 2 ** (m - j)
    ph = {jj: phi(j, jj, p, I) for jj in range(0, d + 1)}
    tot = 0
    for comp in itertools.product(range(d + 1), repeat=nb):
        if sum(comp) != k: continue
        for c0 in range(nb):
            prod = 1
            for c in range(nb):
                if c != c0: prod *= comb(d, comp[c])
            tot += ph[comp[c0]] * prod
    return tot, ph

if __name__ == '__main__':
    m, k = 4, 8
    for p in (1000000009, 998244353):
        print(f"p = {p} (p mod 4 = {p % 4})")
        for j in (1, 2, 3):
            val, ph = f_block(m, k, j, p)
            print(f"  m={m} k={k} rank-{j} block: f >= {val}   phi_{2**j} = {ph}")
            sys.stdout.flush()
    # sanity: block 3 hook k = 4 (80, 144, 224 per b56-methods)
    for j in (1, 2):
        val, ph = f_block(3, 4, j, 1000000009)
        print(f"  sanity m=3 k=4 rank-{j}: f >= {val}")
    # cap formula k(n-k)C(n,k)/(2(n-1))
    n = 16; print("cap formula m=4 k=8:", k * (n - k) * comb(n, k) // (2 * (n - 1)))
