#!/usr/bin/env python3
"""ctr-search / m9.py: Problem D at m = 9 has the single counting profile (7^9).  NC_7 then says every PAIR of the
nine dim-7 port sets covers a co-flat; a pair holds <= 72 ports, so the co-flat has size 53, 56 or 70.
Build the graph on distinct dim-7 port sets (edges = covering pairs) and find its largest clique."""
import sys, os, json, time, pickle, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
from nd import rref_subspaces, mask
t0 = time.time()
cof = [c for c in pickle.load(open(os.path.join(HERE, 'out/coflats80.pkl'), 'rb')) if bin(c).count('1') <= 72]
S7 = set()
for rows in rref_subspaces(9, 2):
    m = (1 << N) - 1
    for r in rows: m &= mask[r]
    S7.add(m)
S7 = [s for s in S7 if bin(s).count('1') >= 53 - 36]
print('dim-7 port sets:', len(S7), 'co-flats <= 72:', len(cof), flush=True)
edges = {}
for c in cof:
    sz = bin(c).count('1')
    big = [k for k, s in enumerate(S7) if bin(s & c).count('1') >= sz - 36]
    for x, y in itertools.combinations(big, 2):
        if c & ~(S7[x] | S7[y]) == 0:
            edges.setdefault(x, set()).add(y); edges.setdefault(y, set()).add(x)
ne = sum(len(v) for v in edges.values()) // 2
print('covering pairs (edges):', ne, f'({time.time()-t0:.0f} s)', flush=True)
# max clique by simple branch and bound
best = [0]
def bk(R, P):
    if not P:
        best[0] = max(best[0], len(R)); return
    if len(R) + len(P) <= best[0]: return
    for v in sorted(P):
        bk(R + [v], P & edges.get(v, set()))
        P = P - {v}
        if len(R) + len(P) <= best[0]: return
bk([], set(edges.keys()))
print('largest clique of dim-7 sets with every pair covering a co-flat:', best[0], '(need 9)', f'({time.time()-t0:.0f} s)')
json.dump(dict(sets=len(S7), edges=ne, max_clique=best[0]), open(os.path.join(HERE, 'm9.json'), 'w'))
