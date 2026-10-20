#!/usr/bin/env python3
"""placements.py (w8-rmexact): classes of cyclic address groups H of GL(4,2) RELATIVE TO THE TARGET I, i.e. orbits of
the orthogonal group O_4 (stabiliser of I under S -> g S g^T) acting on subgroups by conjugation.  Two subgroups in
one class give the same no-phase game up to relabelling (conjugation by g in O_4 fixes the target I).
Matrices: lists of row bitmasks (as forms.py).  Prints, per order, the class sizes and one generator per class."""
import itertools, sys
m = 4; I = [1 << i for i in range(m)]
def mm(A, B):
    out = []
    for i in range(m):
        r = 0
        for k in range(m):
            if (A[i] >> k) & 1: r ^= B[k]
        out.append(r)
    return out
def tr(A): return [sum(((A[r] >> c) & 1) << r for r in range(m)) for c in range(m)]
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
GL = [list(R) for R in itertools.product(range(16), repeat=m) if rank(R) == m]
assert len(GL) == 20160
key = lambda A: tuple(A)
O4 = [g for g in GL if mm(g, tr(g)) == I]
inv = {key(g): next(h for h in GL if mm(g, h) == I) for g in O4}
def order(g):
    h = g; k = 1
    while h != I: h = mm(g, h); k += 1
    return k
def cyc(g):
    grp = [I]; h = g
    while h != I: grp.append(h); h = mm(g, h)
    return frozenset(key(x) for x in grp)
def fixed(g): return any(sum((bin(r & v).count('1') & 1) << i for i, r in enumerate(g)) == v for v in range(1, 16))
want = [int(v) for v in sys.argv[1].split(',')] if len(sys.argv) > 1 else [3, 5, 15]
for n in want:
    subs = {}
    for g in GL:
        if order(g) == n: subs.setdefault(cyc(g), g)
    left = set(subs); classes = []
    while left:
        s0 = next(iter(sorted(left)))
        orb = set()
        for o in O4:
            oi = inv[key(o)]
            orb.add(frozenset(key(mm(mm(o, list(h)), oi)) for h in s0))
        classes.append((len(orb), subs[s0], fixed(subs[s0])))
        left -= orb
    print('order', n, 'subgroups', len(subs), 'classes', len(classes))
    for sz, g, fx in sorted(classes, key=lambda c: (c[2], -c[0])):
        print('   class size %3d  gen rows %s  %s' % (sz, g, 'has fixed vectors' if fx else 'fixed-point-free'))
