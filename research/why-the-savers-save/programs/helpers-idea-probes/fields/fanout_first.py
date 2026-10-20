#!/usr/bin/env python3
"""g-fields I1 relaxed: re-timing allowed. A target's frame when it reads packet P contains every plane it has read before.
For a copy to pass from consumer a to a later consumer b of P, every plane a has read up to P must lie in u_b^perp.
Count, for consumer pairs of the same packet, the planes of a (other than P) that lie in u_b^perp: if none, a must
read P FIRST; each target has only one first read, so the copies one-chain fan-out can save are bounded by a matching.
Usage: python3 -I fanout_first.py"""
import json
from collections import Counter, defaultdict
T = json.load(open('timeline.json')); R = json.load(open('consumers.json'))
ev = {int(q): E for q, E in T['ev'].items()}; cls = {int(q): k for q, k in T['cls'].items()}
ports = T['ports']; v = 112
dot = lambda a, b: bin(a & b).count('1') & 1
inperp = lambda P, t: all(not dot(x, ports[t - v]) for x in P)
cons = defaultdict(set); planes_of = defaultdict(set)
for i, r in enumerate(R):
    for q, k, tv in r['visitors']:
        for (t, ph, role, p, f) in ev[q]:
            if t > tv and role == 'r' and cls[p] == 'y' and 11003 <= t < 12921: cons[i].add(p); planes_of[p].add(i)
print('targets consuming planes in B:', len(planes_of), 'planes per target', Counter(len(s) for s in planes_of.values()))
pair = Counter(); must_first = defaultdict(set)
for i, ts in cons.items():
    for a in ts:
        for b in ts:
            if a == b: continue
            others = [j for j in planes_of[a] if j != i]
            k = sum(1 for j in others if inperp(R[j]['plane'], b))
            pair[(len(others), k)] += 1
print('(other planes of a, of which inside u_b^perp) over ordered consumer pairs:', dict(pair))
# bound: a packet with k consumers served by c copies: each copy chain of length L needs L-1 passers, each passer reads
# P with frame inside the next consumer's perp; with no other plane inside, the passer must read P first.
