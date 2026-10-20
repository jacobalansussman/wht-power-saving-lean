# v-ex2-duality wd_cert.py: is the spreading relation Rs = 2[2] <= [1] + [3] excluded in the WIDENED game Wd
# (free letters: every address permutation, every Z4 phase pattern, Paulis; any mixing)?  Certificate A of w-bounds
# needs Clifford free letters (it permutes Pauli labels), so it is not valid in Wd.  Here: the Wd-valid certificates
# V^(x)k with Y_0 = span{e_x^(x)k} (kept by every monomial matrix, hence by every Wd free letter; g acts as g^(x)k, a
# representation, so Theorem R applies).  phi_r = f(rank-r block at width r); at width w, f = 2^(w-r) phi_r.
# Exact ranks over F_p (p = 1 mod 4, i = sqrt(-1) mod p; matrices scaled by 2 per direction) for two primes,
# plus a float SVD rank.  rank mod p <= true rank (so mod-p values are LOWER bounds); SVD corroborates.
import sys
import numpy as np
from itertools import product
from functools import reduce

def rank_mod_p(rows, p):
    M = [list(r) for r in rows]; rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], p - 2, p)
        M[rk] = [(v * inv) % p for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk

def dir_int(p, isq):  # 2*dir = (1+i) I + (1-i) X  over F_p
    a, b = (1 + isq) % p, (1 - isq) % p
    return [[a, b], [b, a]]
def kron_int(A, B, p):
    n, m = len(A), len(B)
    return [[(A[i // m][j // m] * B[i % m][j % m]) % p for j in range(n * m)] for i in range(n * m)]
def matvec_cols(g):  # columns g e_x
    n = len(g); return [[g[r][x] for r in range(n)] for x in range(n)]

def phi(r, k, p, isq):
    g = reduce(lambda A, B: kron_int(A, B, p), [dir_int(p, isq)] * r)   # 2^r * (rank-r block on r rows)
    n = 2 ** r
    cols = matvec_cols(g)
    # Y_0 basis: e_x^(x)k (index x*(n^(k-1)+...+1)); g Y_0 basis: (g e_x)^(x)k.  dim(Y0 + gY0) - dim Y0.
    def tens(v):
        out = v
        for _ in range(k - 1): out = [(a * b) % p for a in out for b in v]
        return out
    N = n ** k
    Y0 = []
    for x in range(n):
        e = [0] * n; e[x] = 1; Y0.append(tens(e))
    gY0 = [tens(c) for c in cols]
    return rank_mod_p(Y0 + gY0, p) - rank_mod_p(Y0, p)

def phi_float(r, k):
    D = np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]]) / 2
    g = reduce(np.kron, [D] * r); n = 2 ** r
    def tens(v):
        out = v
        for _ in range(k - 1): out = np.kron(out, v)
        return out
    Y0 = [tens(np.eye(n)[x]) for x in range(n)]; gY0 = [tens(g[:, x]) for x in range(n)]
    s = np.linalg.svd(np.array(Y0 + gY0), compute_uv=False)
    return int((s > 1e-9 * s[0]).sum()) - n

for k in (2, 3):
    for r in range(1, 5 if k == 2 else 4):
        vals = [phi(r, k, p, isq) for p, isq in ((1000000009, None), (998244353, None)) for isq in [pow(  # sqrt(-1)
            next(a for a in range(2, 50) if pow(a, (p - 1) // 2, p) == p - 1), (p - 1) // 4, p)]]
        print(f"k={k}: phi_{r} (rank-{r} block at width {r}) mod two primes {vals}, float {phi_float(r, k)}", flush=True)
