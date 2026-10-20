#!/usr/bin/python3
"""even.py (g-engine): the even-weight certificate on 4 columns at block size 2, end to end, exact in F_p.
Y = span of e_x1 (x) e_x2 (x) e_x3 (x) e_x4 with x1+x2+x3+x4 = 0 (addresses of 2 bits).  dim 64 of 256."""
import sys, itertools, random
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from check import *
m = 2; n = 4; N = n ** 4
def idx(x): return ((x[0] * n + x[1]) * n + x[2]) * n + x[3]
Yset = [x for x in itertools.product(range(n), repeat=4) if x[0] ^ x[1] ^ x[2] ^ x[3] == 0]
inY = set(idx(x) for x in Yset)
def apply4(g, x):
    """g^(x4) e_x as a sparse dict"""
    out = {0: 1}
    for c in range(4):
        new = {}
        for k, v in out.items():
            for y in range(n):
                a = g[y][x[c]]
                if a: new[k * n + y] = (new.get(k * n + y, 0) + v * a) % P
        out = new
    return out
def F4(g):
    rows = []
    for x in Yset:
        d = apply4(g, x); r = [0] * N
        for k, v in d.items(): r[k] = v
        # quotient by Y: drop the coordinates in Y
        rows.append([r[j] for j in range(N) if j not in inY])
    return rank(rows)
def keeps(g):
    return all(all((k in inY) or v == 0 for k, v in apply4(g, x).items()) for x in Yset)
print("dim Y =", len(Yset), "of", N)
for zs, name in (((1,), "block rank 1, direction e1"), ((3,), "block rank 1, direction e1+e2"), ((1, 2), "kernel (rank 2)")):
    print(f"  {name}: dim(Y + g^4 Y) - dim Y = {F4(block(m, list(zs)))}; predicted {64 - 64 // 2 ** len(zs)}")
def perm_lin(cols):
    A = zeros(n, n)
    for x in range(n):
        y = 0
        for i in range(m):
            if x >> i & 1: y ^= cols[i]
        A[x][y] = 1
    return A
def phase(z, c): return [[(I if (bin(z & x).count("1") + c) & 1 else 1) if x == y else 0 for y in range(n)] for x in range(n)]
frees = [shift(m, z) for z in (1, 2, 3)] + [perm_lin(c) for c in ((1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2))] + [phase(z, c) for z in range(4) for c in (0, 1)]
print("  all", len(frees), "free letters of block size 2 keep Y:", all(keeps(b) for b in frees))
# column permutations keep Y trivially (the condition x1+x2+x3+x4 = 0 is symmetric); a transposition of ADDRESSES that is not affine does not:
t = eye(n); t[0], t[1] = t[1], t[0]
print("  CONTROL: a non-affine address permutation (swap addresses 0 and 1 only) keeps Y:", keeps(t), "(expected False)")
print("  CONTROL: dir(e1) keeps Y:", keeps(dirm(m, 1)), "(expected False)")
