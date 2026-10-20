#!/usr/bin/env python3
"""s8_hyperplanes.py (t2-block3, stdlib): PROOF.md S3' says a symmetry element can fix a move only if its linear part
on L_Z fixes a hyperplane.  Count, for every element of GL(m,2) (m = 3, 4), by order and conjugacy-invariant data
(order, number of fixed vectors), the number of invariant hyperplanes.  Elements with none have only free move orbits.
"""
import itertools
from collections import defaultdict

def mats(m):
    for bits in itertools.product([0, 1], repeat=m * m):
        A = [bits[i * m:(i + 1) * m] for i in range(m)]
        yield A

def apply(A, x, m):
    y = 0
    for i in range(m):
        s = 0
        for j in range(m):
            s ^= A[i][j] & ((x >> j) & 1)
        y |= s << i
    return y

for m in (3, 4):
    N = 1 << m
    table = defaultdict(lambda: [0, set()])
    for A in mats(m):
        img = [apply(A, x, m) for x in range(N)]
        if len(set(img)) != N:
            continue
        # order
        k, cur = 1, img[:]
        while cur != list(range(N)):
            cur = [img[c] for c in cur]; k += 1
        fixed_vec = sum(1 for x in range(1, N) if img[x] == x)
        # hyperplanes = kernels of nonzero functionals f (f.x = popcount(f & x) mod 2); H_f invariant iff f o A = f on... iff A H_f = H_f
        inv_h = 0
        for f in range(1, N):
            H = [x for x in range(N) if bin(f & x).count("1") % 2 == 0]
            if set(img[x] for x in H) == set(H):
                inv_h += 1
        key = (k, fixed_vec)
        table[key][0] += 1
        table[key][1].add(inv_h)
    print("GL(%d,2):" % m)
    for (k, fv), (cnt, hs) in sorted(table.items()):
        print("   order %2d, fixed nonzero vectors %2d: %5d elements, invariant hyperplanes %s%s" % (
            k, fv, cnt, sorted(hs), "   <- no fixed move possible" if hs == {0} else ""))
