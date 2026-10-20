#!/usr/bin/env python3
"""v_exact_local.py (v-e3-slice6): exact (rational) ranks for the local-problem witnesses claimed by e3-slice6:
 bent circulant on AG(4,2) (L(4,4) feasible), the n=5 circulant from out_circ5.txt (L(5,8) feasible), and the
 Singer permutation on AG(3,2) (L(3,2) feasible). Plus: the plain 1 in 5 satisfies the block-5 star rows."""
from fractions import Fraction as Fr
import itertools
def qrank(M):
    M = [[Fr(x) for x in r] for r in M]; r = 0; rows = len(M); cols = len(M[0]) if rows else 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r
pc = lambda x: bin(x).count('1')
def hyperplanes(n):
    out = []
    for l in range(1, 1 << n):
        for c in (0, 1): out.append([x for x in range(1 << n) if pc(x & l) % 2 == c])
    return out
def circ_from_spectrum(n, s):   # E_xy = sum_u s(u) (-1)^{u.(x+y)}  (= 2^n times the normalised circulant)
    N = 1 << n
    c = [sum(s[u] * (-1) ** pc(u & v) for u in range(N)) for v in range(N)]
    return [[c[x ^ y] for y in range(N)] for x in range(N)]
# bent: E_xy = (-1)^phi(x+y), phi = x1x2 + x3x4
phi = lambda x: ((x & 1) & ((x >> 1) & 1)) ^ (((x >> 2) & 1) & ((x >> 3) & 1))
Eb = [[(-1) ** phi(x ^ y) for y in range(16)] for x in range(16)]
hr = [qrank([[Eb[x][y] for y in H] for x in H]) for H in hyperplanes(4)]
print("bent n=4: rank", qrank(Eb), "; hyperplane ranks min/max", min(hr), max(hr), "over", len(hr))
# n=5 circulant from the truth table printed in out_circ5.txt (index u = position; convention irrelevant up to GL(5,2))
tt = "01101111100100000011101011000101"
s = [(-1) ** int(ch) for ch in tt]
E5 = circ_from_spectrum(5, s)
hr5 = [qrank([[E5[x][y] for y in H] for x in H]) for H in hyperplanes(5)]
ac = [sum(s[u] * s[u ^ h] for u in range(32)) for h in range(1, 32)]
print("circ n=5: rank", qrank(E5), "; hyperplane ranks min/max", min(hr5), max(hr5), "over", len(hr5), "; max nonzero autocorr", max(ac))
# also try the bit-reversed reading of the truth table
tt2 = [tt[int(format(u, '05b')[::-1], 2)] for u in range(32)]
E5b = circ_from_spectrum(5, [(-1) ** int(ch) for ch in tt2])
print("circ n=5 (bit-reversed reading): rank", qrank(E5b), "; max hyperplane rank", max(qrank([[E5b[x][y] for y in H] for x in H]) for H in hyperplanes(5)))
# Singer permutation on F_8 = F_2[a]/(a^3+a+1): x -> a x (nonzero), 0 -> 0 ; E = permutation matrix
def mul_a(x):
    x <<= 1
    if x & 8: x ^= 0b1011
    return x
Ps = [[1 if (y == mul_a(x)) else 0 for y in range(8)] for x in range(8)]
hr3 = [qrank([[Ps[x][y] for y in H] for x in H]) for H in hyperplanes(3)]
print("Singer permutation n=3: rank", qrank(Ps), "; hyperplane ranks min/max", min(hr3), max(hr3))
# plain 1 in 5 vs the block-5 star rows (all lines of all even 3-flat supports), with w_P set to its largest allowed value
m = 5; pts = list(range(1, 32)); mons = [(j, k) for j in range(m) for k in range(j, m)]
monmask = [sum(1 << z for z in pts if (z >> j) & 1 and (z >> k) & 1) for (j, k) in mons]
odd8, even8 = set(), set(); mask = 0; tr = 0
for g in range(1 << len(mons)):
    if g:
        b = (g & -g).bit_length() - 1; mask ^= monmask[b]; tr ^= (mons[b][0] == mons[b][1])
    if pc(mask) == 8: (odd8 if tr else even8).add(mask)
c = {z: (1 if pc(z) == 1 else 0) for z in pts}
def wbound(P):
    Pm = sum(1 << z for z in P)
    return min([sum(c[z] for z in pts if (G >> z) & 1) - 1 for G in odd8 if G & Pm == Pm] + [4 * 1])  # rank/W <= |P| moves/W
worst = None; nrows = 0
for F in even8:
    Fp = [z for z in pts if (F >> z) & 1]
    for x, y in itertools.combinations(Fp, 2):
        planes = []
        for w in Fp:
            if w in (x, y): continue
            P = frozenset((x, y, w, x ^ y ^ w))
            if P <= set(Fp) and P not in planes: planes.append(P)
        assert len(planes) == 3
        slack = sum(min(wbound(P), sum(c[z] for z in P)) for P in planes) - 2 * (c[x] + c[y])
        nrows += 1; worst = slack if worst is None else min(worst, slack)
print("plain 1 in 5: star rows checked", nrows, "; min slack (>= 0 means feasible)", worst, "; even 3-flat supports", len(even8), "; odd", len(odd8))
