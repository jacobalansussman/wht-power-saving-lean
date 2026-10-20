#!/usr/bin/env python3
"""cx-mirror / vclique.py: independent check (plain nested search, no branch and bound) that no 5 x's of any U8
support are pairwise non-orthogonal, i.e. the natural mirror's longest drop-free run inside a support is 4."""
import json, os
D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'all-twists', 'data.json')))
U = D['ports']; dot = lambda a, b: bin(a & b).count('1') & 1
res = []
for S in D['supports']:
    n = len(S); adj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and dot(U[S[i]], U[S[j]]) == 1: adj[i] |= 1 << j
    c4 = c5 = 0
    for a in range(n):
        A = adj[a] & ~((1 << (a + 1)) - 1)
        while A:
            b = (A & -A).bit_length() - 1; A &= A - 1
            B = A & adj[b]
            while B:
                c = (B & -B).bit_length() - 1; B &= B - 1
                C = B & adj[c]
                while C:
                    d = (C & -C).bit_length() - 1; C &= C - 1
                    c4 += 1
                    c5 += bin(C & adj[d]).count('1')
    res.append((c4, c5))
print('per support: (non-orthogonal 4-sets, 5-sets):', res)
