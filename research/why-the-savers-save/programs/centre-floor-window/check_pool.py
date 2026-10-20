#!/usr/bin/env python3
"""ctr-window / check_pool.py POOL: independent premise check for verify_mu. Every pool set must be 8 distinct P112
lines whose roots form a basis of R^8 (exact rank over Q, integer Bareiss determinant != 0); the labels line must
satisfy u_t.u_s = 1 iff t = s or <r_t,r_s> = 0; the first F sets must be pairwise-orthogonal frames."""
import sys, os
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from designs_cx import design, ip
L = design('P112')
lines = open(sys.argv[1]).read().split('\n')
n, P, F = map(int, lines[0].split()); labs = list(map(int, lines[1].split()))
assert n == 112 and len(labs) == 112 and len(set(labs)) == 112
for i in range(n):
    for j in range(n):
        assert (bin(labs[i] & labs[j]).count('1') & 1) == (1 if (i == j or ip(L[i], L[j]) == 0) else 0)
def det(M):            # Bareiss, exact integers
    M = [row[:] for row in M]; k = len(M); prev = 1; sign = 1
    for c in range(k - 1):
        if M[c][c] == 0:
            sw = next((r for r in range(c + 1, k) if M[r][c] != 0), None)
            if sw is None: return 0
            M[c], M[sw] = M[sw], M[c]; sign = -sign
        for r in range(c + 1, k):
            for s in range(c + 1, k):
                M[r][s] = (M[r][s] * M[c][c] - M[r][c] * M[c][s]) // prev
        prev = M[c][c]
    return sign * M[k - 1][k - 1]
sets = [list(map(int, l.split())) for l in lines[2:] if l.strip()]
assert len(sets) == P
for idx, S in enumerate(sets):
    assert len(S) == 8 and len(set(S)) == 8 and all(0 <= t < n for t in S)
    assert det([list(L[t]) for t in S]) != 0, ('not a basis', idx)
    if idx < F: assert all(ip(L[a], L[b]) == 0 for a in S for b in S if a != b), ('frame not orthogonal', idx)
print('pool OK:', P, 'sets, all bases of R^8;', F, 'of them orthogonal frames; labels verified')
