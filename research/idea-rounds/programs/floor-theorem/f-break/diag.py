#!/usr/bin/env python3
"""diag.py: exact check (modulo two primes) of the transfer ranks used by the twin-column floor.
V = C^d (x) C^d, d = 2^m, two digit columns.  D = span{ e_x (x) e_x }.  For a matrix g on one column,
q(g) = rank of  D -> V/D,  v |-> (g (x) g) v mod D.   Claim: q(block of rank k) = 2^m - 2^(m-k); q(kernel) = 2^m - 1;
q(monomial) = 0.  Also checks a random monomial map and a shift on ONE column only (which must NOT preserve D)."""
import sys, random
def rank(rows, p):
    rows = [r[:] for r in rows]; rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = pow(rows[rk][c], p-2, p); rows[rk] = [x*iv % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]; rows[i] = [(x - f*y) % p for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk
def kron(A, B, p):
    return [[A[i][j]*B[k][l] % p for j in range(len(A)) for l in range(len(B))] for i in range(len(A)) for k in range(len(B))]
def q(g1, g2, d, p):
    # images of e_x (x) e_x under g1 (x) g2, with the diagonal coordinates (y,y) deleted
    out = []
    for x in range(d):
        out.append([g1[y1][x]*g2[y2][x] % p for y1 in range(d) for y2 in range(d) if y1 != y2])
    return rank(out, p)
for p in (65537, 40961):
    I = next(x for x in range(2, p) if x*x % p == p-1); h = pow(2, p-2, p)
    C = [[(1+I)*h % p, (1-I)*h % p], [(1-I)*h % p, (1+I)*h % p]]; I2 = [[1, 0], [0, 1]]
    for m in (1, 2, 3, 4, 5):
        d = 2**m; res = []
        for k in range(0, m+1):
            g = [[1]]
            for b in range(m): g = kron(g, C if b < k else I2, p)
            res.append(q(g, g, d, p))
        want = [d - d//2**k for k in range(m+1)]
        random.seed(m); perm = list(range(d)); random.shuffle(perm)
        mono = [[0]*d for _ in range(d)]
        for x in range(d): mono[perm[x]][x] = random.randrange(1, p)
        ident = [[int(i == j) for j in range(d)] for i in range(d)]
        shift = [[int(i == j ^ 1) for j in range(d)] for i in range(d)]
        print(f"p={p} m={m}: q(rank k block), k=0..m = {res}  claimed {want}  match={res == want};"
              f"  q(random monomial on both columns)={q(mono, mono, d, p)};  q(shift on column 1 only)={q(shift, ident, d, p)}", flush=True)
