#!/usr/bin/env python3
"""v_nearU6.py (v-e3-slice6): my own recount of the near-uniform exclusions at block 6 (U_6 plus a multiset Z of t
extra directions, W arrays). Criterion (Theorem C applied to one even support): some EVEN weight-16 support F with
no extra move in it, such that every one of F's 30 affine hyperplanes H lies in an ODD support G with
|D_G| - W = 16 + #(Z in G, with multiplicity) - W <= 3. (Weight-16 supports are the only ones that can have slack <= 3.)
"""
import numpy as np, itertools, sys
m = 6; pts = list(range(1, 64))
mons = [(j, k) for j in range(m) for k in range(j, m)]
monmask = [sum(1 << z for z in pts if (z >> j) & 1 and (z >> k) & 1) for (j, k) in mons]
odd, even = [], []
mask = 0; tr = 0
for g in range(1 << len(mons)):
    if g:
        b = (g & -g).bit_length() - 1; mask ^= monmask[b]; tr ^= (mons[b][0] == mons[b][1])
    if bin(mask).count('1') == 16: (odd if tr else even).append(mask)
odd = sorted(set(odd)); even = sorted(set(even))
mem = lambda S: np.array([[(s >> z) & 1 for z in range(64)] for s in S], dtype=np.int64)
Gm = mem(odd); Fm = mem(even)
# hyperplanes of each even F and their odd covers
cov = np.full((len(even), 30, 7), len(odd), dtype=np.int64)
for fi, F in enumerate(even):
    Fp = [z for z in pts if (F >> z) & 1]; a = Fp[0]
    hy = set()
    for l in range(1, 64):
        for c in (0, 1):
            H = sum(1 << z for z in Fp if bin(z & l).count('1') % 2 == c)
            if bin(H).count('1') == 8: hy.add(H)
    assert len(hy) == 30, len(hy)
    for hi, H in enumerate(sorted(hy)):
        gs = [gi for gi, G in enumerate(odd) if G & H == H]
        assert 0 <= len(gs) <= 7, len(gs)
        cov[fi, hi, :len(gs)] = gs
def excluded(zvec, W):
    cnt = Gm @ zvec
    ok = np.append(16 + cnt - W <= 3, False)
    good = ok[cov].any(axis=2).all(axis=1) & ((Fm @ zvec) == 0)
    return bool(good.any())
for W, t in [(13, 0), (13, 1), (13, 2), (14, 1), (14, 2), (14, 3), (15, 1), (15, 2), (15, 3), (16, 1), (16, 2), (16, 3)]:
    n = ex = 0; openl = []
    for Z in itertools.combinations_with_replacement(pts, t):
        zvec = np.zeros(64, dtype=np.int64)
        for z in Z: zvec[z] += 1
        n += 1
        if excluded(zvec, W): ex += 1
        elif len(openl) < 200: openl.append(Z)
    allfull = all(63 in Z or any(a ^ b == 63 for a, b in itertools.combinations(Z, 2)) for Z in openl)
    print("W=%d t=%d: multisets %d, excluded %d, open %d%s" % (W, t, n, ex, n - ex,
          ("; every open case contains 63 or two extras XOR to 63: %s; e.g. %s" % (allfull, openl[:4])) if n > ex else ""))
    sys.stdout.flush()
