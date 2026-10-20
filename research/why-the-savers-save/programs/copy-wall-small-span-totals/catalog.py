"""catalog.py (d-R7): sub-units of P112 = ports whose labels lie in one hyperplane (dim 8) or one 7-dim subspace of
F_2^9.  For each distinct port set: size, rank of roots, largest orthogonal clique, count of hyperplanes giving it."""
import sys, time, itertools, collections; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import warnings; warnings.filterwarnings('ignore')
from model import *
import numpy as np
t0 = time.time()
OA = [sum(1 << j for j in range(NP) if ORTH[i, j]) for i in range(NP)]
def omega(us):
    best = 0
    def rec(size, cand):
        nonlocal best
        if size > best: best = size
        c = cand
        while c:
            if size + bin(c).count('1') <= best: return
            low = c & -c; i = low.bit_length() - 1; c ^= low
            rec(size + 1, c & OA[i])
    rec(0, sum(1 << i for i in us)); return best
seen = collections.Counter(); info = {}
for h in range(1, 512):
    us = tuple(lines_of_dual([h]))
    seen[us] += 1
    if us not in info: info[us] = (len(us), qrank(W[list(us)]) if us else 0, omega(us), h)
print('hyperplanes: distinct port sets', len(seen))
for (n, r, om, h), c in sorted(collections.Counter(info[u] [:3] for u in seen for _ in range(seen[u])).items(), reverse=True) if False else []: pass
tab = collections.Counter((info[u][0], info[u][1], info[u][2]) for u in seen for _ in range(seen[u]))
for k in sorted(tab, reverse=True): print('  hyperplane sub-unit (ports, root rank, omega) %s : %d hyperplanes' % (k, tab[k]))
# example h per type
ex = {}
for u in seen:
    k = info[u][:3]; ex.setdefault(k, info[u][3])
print('  example normals h per type:', {k: ex[k] for k in sorted(ex, reverse=True)[:8]})
# 7-dim subspaces with the most ports
best = collections.Counter(); exs = {}
for h1 in range(1, 512):
    for h2 in range(h1 + 1, 512):
        if (h1 ^ h2) < h2: continue
        us = tuple(lines_of_dual([h1, h2]))
        if len(us) >= 30:
            k = (len(us), qrank(W[list(us)]), omega(us)); best[k] += 1; exs.setdefault(k, (h1, h2))
for k in sorted(best, reverse=True): print('  7-dim sub-unit (ports, rank, omega) %s : %d subspaces, e.g. duals %s' % (k, best[k], exs[k]))
print('%.0fs' % (time.time() - t0))
