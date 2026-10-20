#!/usr/bin/env python3
"""g-fields I6 (circuit CSE / shared intermediate frames): relax H2 (exact output planes). A packet row may be created at
any frame F with P <= F <= S_P, where S_P = intersection of u_t^perp over P's B-phase consumers (frame sandwich).
Then ONE register can be a creation visitor for two planes P1 then P2 if P1 + P2 <= S_P2 (its chain P1 <= F1 <= F2,
F2 >= P1 + P2). Visit toy: each plane needs 3 visits (two copies), sources give their one visit as before; a helper gives
one visit, or two visits to a compatible ordered pair. Least helpers = 392 - (max number of disjoint pair-slots), bounded
here by a greedy b-matching (lower estimate of the saving) and by the slot count (upper). Necessary condition only.
Usage: python3 -I relaxH2.py"""
import json
from collections import Counter, defaultdict
T = json.load(open('timeline.json')); R = json.load(open('consumers.json'))
ev = {int(q): E for q, E in T['ev'].items()}; cls = {int(q): k for q, k in T['cls'].items()}
ports = T['ports']; v = 112
dot = lambda a, b: bin(a & b).count('1') & 1
cons = defaultdict(set)
for i, r in enumerate(R):
    for q, k, tv in r['visitors']:
        for (t, ph, role, p, f) in ev[q]:
            if t > tv and role == 'r' and cls[p] == 'y' and 11003 <= t < 12921: cons[i].add(p)
n = len(R)
inS = lambda P, j: all(not dot(x, ports[t - v]) for t in cons[j] for x in P)
helper_slots = {i: 3 - sum(1 for q, k, tv in r['visitors'] if k == 'x') for i, r in enumerate(R)}
print('helper visits needed', sum(helper_slots.values()))
E = [(a, b) for a in range(n) for b in range(n) if a != b and inS(R[a]['plane'], b)]
und = defaultdict(set)
for a, b in E: und[a].add(b); und[b].add(a)
print('ordered compatible pairs', len(E), 'planes with >= 1 partner', sum(1 for i in range(n) if und[i]),
      'unconsumed planes (S = all)', sum(1 for i in range(n) if not cons[i]))
# exact b-matching via max-flow on the undirected graph (each plane offers helper_slots[i] visit-slots; a pair uses one slot at each end)
# general-graph b-matching is not bipartite; use the bipartite double cover bound (upper bound) and a greedy (lower).
slots = dict(helper_slots); used = 0
for a, b in sorted({tuple(sorted(e)) for e in E}, key=lambda e: len(und[e[0]]) + len(und[e[1]])):
    if slots[a] > 0 and slots[b] > 0: slots[a] -= 1; slots[b] -= 1; used += 1
print('greedy disjoint pair-slots:', used, '-> helpers', sum(helper_slots.values()) - used, '(from 392)')
print('upper bound on pairs: floor(slots/2) =', sum(helper_slots.values()) // 2)
# with consumers-only sandwich, also note how many pairs are within one block vs across blocks
lab = lambda i: frozenset(s for s in range(v) if any(ports[s] == x for x in []))
print('pair consumer-count profile:', Counter((len(cons[a]), len(cons[b])) for a, b in E).most_common(8))
