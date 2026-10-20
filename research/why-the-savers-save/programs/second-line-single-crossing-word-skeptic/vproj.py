#!/usr/bin/env python3
"""v-an2-word: the projection identity, own code.  For an elementary paid step (a climb f < g or a drop g < f) between
subspaces of F_2^m and the shared line s, cost(f -> g) - cost(proj f -> proj g) = [s in exactly one of f, g], where
proj f = (f + <s>)/<s>.  Algebra: dim proj f = dim f - [s in f], so for nested f, g the difference is exactly
|[s in g] - [s in f]|.  Random check with frames that contain s + a vectors (not just coordinate frames), m = 9, and a
check that the projection keeps equal frames equal and maps the two-sided endpoint frames to the one-sided ones.
usage: python3 -I vproj.py"""
import random
def basis(rows):
    b = {}
    for r in rows:
        while r:
            h = r.bit_length() - 1
            if h in b: r ^= b[h]
            else: b[h] = r; break
    return b
def dim(rows): return len(basis(rows))
def has(rows, v): return dim(rows + [v]) == dim(rows)
m = 9; S = 1 << (m - 1)            # s = the new (last) coordinate
rnd = random.Random(2026)
bad = cr = 0; N = 60000
for _ in range(N):
    f = [rnd.randrange(1 << m) for _ in range(rnd.randrange(0, m))]
    if rnd.random() < 0.5: f.append(S ^ rnd.randrange(1 << (m - 1)))     # an s + a vector
    if rnd.random() < 0.3: f.append(S)
    g = f + [rnd.randrange(1 << m) for _ in range(rnd.randrange(1, 4))]
    for a, b in ((f, g), (g, f)):
        c = abs(dim(b) - dim(a))
        pc = abs((dim(b + [S]) - 1) - (dim(a + [S]) - 1))
        x = has(a, S) != has(b, S); cr += x
        if c - pc != x: bad += 1
print('elementary steps checked %d, crossings %d, identity failures %d' % (2 * N, cr, bad))
# endpoints: u in old coordinates; X: <u> -> s^perp ; Y: <s> -> u^perp (in F_2^m).  Projections: <u> -> all of F_2^(m-1);
# 0 -> u^perp in F_2^(m-1).
u = 0b10110101 & ((1 << (m - 1)) - 1)
sperp = [1 << i for i in range(m - 1)]
uperp = [v for v in range(1 << m) if bin(v & u).count('1') % 2 == 0]
print('X end s^perp: dim %d, proj dim %d (want %d)' % (dim(sperp), dim(sperp + [S]) - 1, m - 1))
print('Y start <s>: proj dim %d (want 0); Y end u^perp: dim %d contains s %s, proj dim %d (want %d)' % (
      dim([S, S]) - 1, dim(uperp), has(uperp, S), dim(uperp + [S]) - 1, m - 2))
