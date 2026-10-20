#!/usr/bin/env python3
"""tclasses.py (t2-block4): classes of target multisets {I, T2, ..., TW} for the multi-target slice at block m, up to the
orthogonal group (g g^T = I), which fixes I and acts by T -> g T g^T.  T ranges over h h^T (h in GL(m,2)).
Prints, for W = 2, one representative T2 per class with class size; for W = 3, the number of classes."""
import itertools, sys
m = int(sys.argv[1]) if len(sys.argv) > 1 else 4
def par(x): return bin(x).count('1') & 1
def mm(A, B):
    out = []
    for i in range(m):
        r = 0
        for k in range(m):
            if (A[i] >> k) & 1: r ^= B[k]
        out.append(r)
    return tuple(out)
def tr(A): return tuple(sum(((A[r] >> c) & 1) << r for r in range(m)) for c in range(m))
def rank(A):
    rows = list(A); r = 0
    for c in range(m):
        p = next((i for i in range(r, m) if (rows[i] >> c) & 1), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(m):
            if i != r and (rows[i] >> c) & 1: rows[i] ^= rows[r]
        r += 1
    return r
I = tuple(1 << i for i in range(m))
GL = [g for g in itertools.product(range(1, 1 << m), repeat=m) if rank(g) == m]
T = sorted({mm(g, tr(g)) for g in GL})
O = [g for g in GL if mm(g, tr(g)) == I]
print('m=%d |GL|=%d targets=%d |O|=%d' % (m, len(GL), len(T), len(O)))
seen = set(); cls = []
for t in T:
    if t in seen: continue
    orb = {mm(mm(g, t), tr(g)) for g in O}
    seen |= orb; cls.append((t, len(orb)))
for t, sz in cls:
    # minimal number of rank-one terms is m; note whether t is I
    print('T2 rows %s class size %d %s' % ('.'.join(map(str, t)), sz, '(I)' if t == I else ''))
print('pair classes:', len(cls))
