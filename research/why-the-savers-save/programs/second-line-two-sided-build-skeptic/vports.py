#!/usr/bin/env python3
"""v-an-build: F_2 rank of the 112 port labels of the U8X word; x-dependent slots of M; scatter reader graph connectivity."""
import gzip, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vreplay import load, forward
d = load(os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz')))
v, R = d['v'], d['R']; n = 2 * v + R
basis = []
for p in d['ports']:
    for b in basis:
        p = min(p, p ^ b)
    if p: basis.append(p)
print('ports', len(d['ports']), 'distinct', len(set(d['ports'])), 'F2 rank', len(basis), 'odd weight all', all(bin(p).count('1') % 2 for p in d['ports']))
M = forward(d)
print('slots with x-dependence in M:', sum(1 for q in range(2 * v, n) if any(k < v for k in M[q])),
      '; not restored:', sum(1 for q in range(2 * v, n) if M[q] != {q: 1}))
# reader graph: y_T ~ y_T' if they read a common centre copy at the scatter
readers = {}
for T, ent in enumerate(d['scat']['table']):
    for k, a, b in ent: readers.setdefault(k, set()).add(T)
print('centre k -> number of y readers:', {k: len(s) for k, s in sorted(readers.items())})
seen = {0}; stack = [0]
while stack:
    T = stack.pop()
    for k, s in readers.items():
        if T in s:
            for U in s:
                if U not in seen: seen.add(U); stack.append(U)
print('reader graph connected:', len(seen) == v)
