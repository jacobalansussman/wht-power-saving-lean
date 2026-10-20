#!/usr/bin/env python3
"""ctr-search / cover3b.py: profile P_B of Problem D at m = 10 needs >= 7 directions of label dim <= 6 such that EVERY
three of their port sets cover a co-flat; by cover3.py only the co-flats N_i can be covered by three dim-6 sets.
Build the 3-uniform hypergraph of covering triples (all 8 N_i, sets = port sets S_W of dim-6 subspaces, taken as
maximal: S_W itself, not its intersection) and look for 4 sets with all 4 triples covering (necessary for 7).
Dim < 6 subspaces are contained in dim-6 ones, so using dim exactly 6 is WLOG for covering."""
import sys, os, itertools, json, time, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
from nd import rref_subspaces, mask
from cover3 import rank_of_mask
t0 = time.time()
Ni = [sum(1 << b for b in range(N) if LINES[b][i] != 0) for i in range(8)]
# all distinct dim-6 port sets with >= 16 ports inside some N_i
S6 = {}
for rows in rref_subspaces(9, 3):
    m = (1 << N) - 1
    for r in rows: m &= mask[r]
    if any(bin(m & g).count('1') >= 16 for g in Ni): S6[m] = rows
sets = list(S6.keys())
print('dim-6 port sets with >= 16 inside some N_i:', len(sets), f'({time.time()-t0:.0f} s)', flush=True)
idx = {s: k for k, s in enumerate(sets)}
triples = set()
for i, g in enumerate(Ni):
    inside = [k for k, s in enumerate(sets) if bin(s & g).count('1') >= 16]
    for x in range(len(inside)):
        for y in range(x + 1, len(inside)):
            a, b = sets[inside[x]], sets[inside[y]]
            u = (a | b) & g
            if bin(u).count('1') < 36: continue
            rest = g & ~u
            # every dim-6 set containing rest: check members of 'inside' (others hold < 16 of N_i but rest may be small)
            for z in range(y + 1, len(inside)):
                c = sets[inside[z]]
                if rest & ~c == 0: triples.add(tuple(sorted((inside[x], inside[y], inside[z]))))
print('covering triples (all N_i):', len(triples), f'({time.time()-t0:.0f} s)', flush=True)
adj = collections.defaultdict(set)
for t in triples:
    for p in itertools.combinations(t, 2):
        adj[p].add(next(v for v in t if v not in p))
# K4^(3): a,b,c,d with abc, abd, acd, bcd all triples
k4 = 0
for (a, b), cs in adj.items():
    cs = sorted(cs)
    for c, d in itertools.combinations(cs, 2):
        if tuple(sorted((a, c, d))) in triples and tuple(sorted((b, c, d))) in triples: k4 += 1
print('4-sets with all four triples covering:', k4, f'({time.time()-t0:.0f} s)', flush=True)
json.dump(dict(sets=len(sets), triples=len(triples), K4=k4), open(os.path.join(HERE, 'cover3b.json'), 'w'), indent=1)
