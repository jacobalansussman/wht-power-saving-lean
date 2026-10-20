#!/usr/bin/env python3
"""my_orbits.py (v-e3-rmbuild): GL(m,2) orbits on binary quadratic forms Q (Q -> Q o g), own code."""
import itertools, sys
m = int(sys.argv[1]); N = 1 << m
def lin(cols):
    img = []
    for x in range(N):
        y = 0
        for k in range(m):
            if x >> k & 1: y ^= cols[k]
        img.append(y)
    return img
GL = [g for g in (lin(c) for c in itertools.product(range(1, N), repeat=m)) if len(set(g)) == N]
pairs = [(i, j) for i in range(m) for j in range(i, m)]
forms = set()
for bits in range(1 << len(pairs)):
    q = [p for b, p in enumerate(pairs) if bits >> b & 1]
    forms.add(tuple(sum((z >> i & 1) & (z >> j & 1) for i, j in q) & 1 for z in range(N)))
left = set(forms); sizes = []
while left:
    Q = left.pop(); orb = {tuple(Q[g[z]] for z in range(N)) for g in GL}
    left -= orb; sizes.append(len(orb))
print("m=%d |GL|=%d forms=%d orbit sizes %s" % (m, len(GL), len(forms), sorted(sizes)))
