#!/usr/bin/env python3
"""v-centre / k4.py: does the 3-uniform hypergraph of covering triples contain 4 vertices with all 4 triples?
Also reports the largest t such that some t vertices have every triple present, by a small exact search."""
import sys, itertools
T = set()
for l in open(sys.argv[1]):
    a, b, c = map(int, l.split()); T.add((a, b, c))
pair = {}
for a, b, c in T:
    for x, y, z in ((a, b, c), (a, c, b), (b, c, a)):
        pair.setdefault((x, y), set()).add(z)
k4 = 0
for (a, b, c) in T:
    # d must form triples with every pair of {a,b,c}
    common = pair[(a, b)] & pair[(a, c)] & pair[(b, c)]
    common.discard(a); common.discard(b); common.discard(c)
    k4 += len(common)
print('triples', len(T), 'vertices', len({v for t in T for v in t}), 'K4 incidences', k4)
