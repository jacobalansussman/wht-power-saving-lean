# Independent check of the group-theoretic inputs of Lemma C at block 2 (BFS mod p), and of L, n_k at m = 2, 3, 4
# by brute-force enumeration of Lagrangian subspaces of F_2^(2m) (m = 2, 3) and the q-formula (m = 4).
import sys, os, itertools, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import prod
from common import letters, mm
p = int(sys.argv[1]) if len(sys.argv) > 1 else 2000029
def bfs(gens, n):
    key = lambda M: M.tobytes()
    e = np.eye(n, dtype=np.int64); seen = {key(e): e}; frontier = [e]
    while frontier:
        nf = []
        for M in frontier:
            for g in gens:
                X = mm(g, M, p); k = key(X)
                if k not in seen: seen[k] = X; nf.append(X)
        frontier = nf
    return seen
free, Wr = letters(2, p)
B = bfs(free, 4); G = bfs(free + [Wr[1]], 4)
print('m=2 |B| =', len(B), ' |G| =', len(G), ' [G:B] =', len(G) / len(B))
# double cosets B W_k B and |B cap W_k B W_k^-1|
Bl = list(B.values())
for k in range(3):
    Wk = Wr[k]
    # B cap W_k B W_k^-1 : b in B with W_k^-1 b W_k in B
    # inverse mod p by Gauss-Jordan
    def invmod(M):
        n = M.shape[0]; A = np.concatenate([M % p, np.eye(n, dtype=np.int64)], 1)
        for c in range(n):
            r = c + np.nonzero(A[c:, c])[0][0]; A[[c, r]] = A[[r, c]]
            A[c] = A[c] * pow(int(A[c, c]), p - 2, p) % p
            for rr in range(n):
                if rr != c and A[rr, c]: A[rr] = (A[rr] - A[rr, c] * A[c]) % p
        return A[:, n:]
    Wki = invmod(Wk)
    inter = sum(1 for b in Bl if mm(mm(Wki, b, p), Wk, p).tobytes() in B)
    print(f'  k={k}: |B W_k B| = |B|^2/|B cap ..| = {len(B)**2//inter} = {len(B)/inter} |B|;  |B cap W_k^-1 B W_k| = {inter} = |B|/{len(B)/inter}')
print('  sum of double coset sizes / |B| should be L = 15')
# Lagrangians in F_2^(2m), symplectic form <(a,b),(c,d)> = a.d + b.c; L0 = {(a,0)}; distance k = m - dim(L cap L0)
def lagr(m):
    N = 2 * m; vecs = range(1, 1 << N)
    def form(u, v):
        a, b = u & ((1 << m) - 1), u >> m; c, d = v & ((1 << m) - 1), v >> m
        return (bin(a & d).count('1') + bin(b & c).count('1')) & 1
    # enumerate m-dim isotropic subspaces by reduced echelon bases: enumerate spans of m independent vectors (sets)
    subs = set()
    def rec(basis, span):
        if len(basis) == m: subs.add(frozenset(span)); return
        for v in vecs:
            if v in span: continue
            if basis and v < basis[-1]: continue
            if all(form(v, w) == 0 for w in basis) and form(v, v) == 0:
                rec(basis + [v], span | {s ^ v for s in span})
    rec([], {0})
    L0 = frozenset(a for a in range(1 << m))
    cnt = {}
    for S in subs:
        k = m - (len(S & L0).bit_length() - 1)
        cnt[k] = cnt.get(k, 0) + 1
    return len(subs), dict(sorted(cnt.items()))
def gb(m, k): return prod(2**(m - i) - 1 for i in range(k)) // prod(2**(i + 1) - 1 for i in range(k))
for m in (2, 3):
    tot, cnt = lagr(m)
    print(f'm={m}: #Lagrangians (brute force) = {tot}, by distance {cnt}; formula prod(2^i+1) = {prod(2**i + 1 for i in range(1, m + 1))}, n_k formula = {[2**(k*(k+1)//2) * gb(m, k) for k in range(m + 1)]}')
m = 4; print('m=4 formula L =', prod(2**i + 1 for i in range(1, 5)), 'n_k =', [2**(k*(k+1)//2) * gb(4, k) for k in range(5)], 'sum', sum(2**(k*(k+1)//2) * gb(4, k) for k in range(5)))
for m in (2, 3, 4):
    Sp = 2**(m*m) * prod(4**i - 1 for i in range(1, m + 1)); L = prod(2**i + 1 for i in range(1, m + 1))
    print(f'm={m}: |G| = 4^(m+1)|Sp| = {4**(m+1) * Sp}, |B| = |G|/L = {4**(m+1) * Sp // L}')
