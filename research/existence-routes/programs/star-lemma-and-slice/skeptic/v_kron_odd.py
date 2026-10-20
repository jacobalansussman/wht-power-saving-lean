#!/usr/bin/env python3
"""v_kron_odd.py (v-e3-slice6): the local-route floor at ODD n. E = B_(n-1) (x) [[0,1],[1,0]], B_(n-1) the bent circulant
on AG(n-1,2) (n-1 even). Claim (mine): nonsingular, rank <= 2^(n-2) on every affine hyperplane of AG(n,2). Exact check at
n = 3, 5. Also the Singer permutation of F_16 (L(4,4) over every field)."""
from fractions import Fraction as Fr
exec(open('v_exact_local.py').read().split('# bent: E_xy')[0])
def bentmat(k):
    def phi(x): return sum(((x >> (2 * i)) & 1) & ((x >> (2 * i + 1)) & 1) for i in range(k // 2)) & 1
    return [[(-1) ** phi(x ^ y) for y in range(1 << k)] for x in range(1 << k)]
for n in (3, 5):
    B = bentmat(n - 1); Nn = 1 << n
    # point (x, b) -> index 2*x + b ; second factor the swap [[0,1],[1,0]]
    E = [[B[i >> 1][j >> 1] * (1 if (i & 1) != (j & 1) else 0) for j in range(Nn)] for i in range(Nn)]
    hr = [qrank([[E[x][y] for y in H] for x in H]) for H in hyperplanes(n)]
    print("n=%d bent(x)swap: rank %d (of %d); hyperplane ranks min %d max %d; bound 2^(n-2) = %d" % (n, qrank(E), Nn, min(hr), max(hr), 1 << (n - 2)))
def mul_a16(x):
    x <<= 1
    if x & 16: x ^= 0b10011
    return x
P = [[1 if y == mul_a16(x) else 0 for y in range(16)] for x in range(16)]
hr = [qrank([[P[x][y] for y in H] for x in H]) for H in hyperplanes(4)]
print("Singer permutation F_16: rank", qrank(P), "; hyperplane ranks min/max", min(hr), max(hr))
