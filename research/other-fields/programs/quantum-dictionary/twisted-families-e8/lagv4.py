#!/usr/bin/env python3
"""lagv4.py (t-family) -- v = 4, every directed pattern, payload F_p for p = 3, 5, 7.
Budget 3 < v needs a term of cost 2 whose F_p-rank of N[R, D] is >= 3 (sum of F_p ranks of the masked terms must be
>= v = 4, a cost-1 term has F_p-rank 1).  Patterns without such a submatrix have floor 4 at once; the others go to
lagmitm.floor_lt_v.  usage: python3 -I lagv4.py p"""
import sys, os, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from lagfloor import all_patterns, rank2
from lagmitm import floor_lt_v
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '../lib')))
from copyfloor import rank
p = int(sys.argv[1]); v = 4; t = time.time()
subs = [c for k in range(1, v + 1) for c in itertools.combinations(range(v), k)]
cnt = cand = hits = 0; ex = None
for N in all_patterns(v):
    cnt += 1
    need = any(rank2([sum(N[T][S] << i for i, S in enumerate(D)) for T in R]) == 2 and
               rank([[N[T][S] for S in D] for T in R], p) >= 3 for R in subs for D in subs)
    if not need:
        continue
    cand += 1
    ok, c = floor_lt_v(N, p)
    if ok:
        hits += 1; ex = ex or N
print('v=4 payload F_%d: %d patterns, %d pass the rank filter, %d with floor < 4 %s (%.1fs)'
      % (p, cnt, cand, hits, ex or '', time.time() - t))
