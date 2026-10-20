#!/usr/bin/env python3
"""orbits_gl.py: orbits of GL(m,2) on binary quadratic forms Q (Q(0) = 0) on F_2^m, by BFS over generators."""
import sys
m = int(sys.argv[1]) if len(sys.argv) > 1 else 4
N = 1 << m
pairs = [(j, k) for j in range(m) for k in range(j, m)]
def table(mask):
    return tuple(sum(1 for b, (j, k) in enumerate(pairs) if (mask >> b) & 1 and (z >> j) & 1 and (z >> k) & 1) & 1 for z in range(N))
tabs = [table(mk) for mk in range(1 << len(pairs))]
idx = {t: i for i, t in enumerate(tabs)}
# generators of GL(m,2): swap e0,e1; cycle of coordinates; transvection e0 -> e0 + e1
def mat_from_images(img):
    return [sum(img[j] for j in range(m) if (z >> j) & 1) and 0 or 0 for z in range(N)]
def lin(img):
    out = []
    for z in range(N):
        v = 0
        for j in range(m):
            if (z >> j) & 1: v ^= img[j]
        out.append(v)
    return out
gens = [lin([2, 1] + [1 << j for j in range(2, m)]), lin([1 << ((j + 1) % m) for j in range(m)]), lin([3] + [1 << j for j in range(1, m)])]
seen = [False] * len(tabs); sizes = []
for i in range(len(tabs)):
    if seen[i]: continue
    stack = [i]; seen[i] = True; cnt = 0
    while stack:
        a = stack.pop(); cnt += 1
        for g in gens:
            t = tuple(tabs[a][g[z]] for z in range(N)); b = idx[t]
            if not seen[b]: seen[b] = True; stack.append(b)
    sizes.append(cnt)
print("m=%d forms=%d orbits=%d sizes=%s" % (m, len(tabs), len(sizes), sorted(sizes)))
