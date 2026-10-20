#!/usr/bin/env python3
"""lagv4b.py (t-family) -- v = 4, every directed pattern, payload F_p: exact test of floor < 4.
A decoder of cost <= 3 must be one term of cost 2 plus one term of cost 1 [proved on paper: a single term fits only
if every pair knows (N = I, cost 4); a cost-1 term has F_p-rank 1 and a cost-2 term F_p-rank <= 3, while the masked
terms must add up to a rank-4 diagonal].  So: for every cost-2 term t1, look up a cost-1 term t2 with
offdiag(t2) = -offdiag(t1) and diag(t1 + t2) nowhere zero.  usage: python3 -I lagv4b.py p"""
import sys, os, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from lagfloor import all_patterns
from lagmitm import terms
p = int(sys.argv[1]); t = time.time(); cnt = hits = 0; ex = None
for N in all_patterns(4):
    cnt += 1
    best, noff = terms(N, p)
    one = {}
    for (od, dg), c in best.items():
        if c == 1: one.setdefault(od, []).append(dg)
    found = False
    for (od, dg), c in best.items():
        if c != 2: continue
        neg = tuple((-x) % p for x in od)
        for dg2 in one.get(neg, ()):
            if all((a + b) % p for a, b in zip(dg, dg2)):
                found = True; break
        if found: break
    if found:
        hits += 1; ex = ex or N
print('v=4 payload F_%d (cost-2 + cost-1 test): %d patterns, %d with floor < 4 %s (%.1fs)' % (p, cnt, hits, ex or '', time.time() - t))
