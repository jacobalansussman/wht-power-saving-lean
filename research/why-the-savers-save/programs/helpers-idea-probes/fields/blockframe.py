#!/usr/bin/env python3
"""g-fields I5 (sparse factorisation / supernode lens): could packets be held at the rank-4 BLOCK frame Q_b (span of the
block's 8 labels) instead of at their 12 distinct rank-3 planes? Then all gates of a block act at one frame, the
packet-plane visit bound (incomparable planes) no longer applies, and the count becomes a register-pressure question.
Necessary condition: every B-phase consumer t of a block-b packet has Q_b inside u_t^perp (frame sandwich), and reads
the packet while its own frame is inside Q_b, so all of t's packets must come from that one block.
Usage: python3 -I blockframe.py"""
import json, itertools
from collections import Counter, defaultdict
T = json.load(open('timeline.json')); R = json.load(open('consumers.json'))
ev = {int(q): E for q, E in T['ev'].items()}; cls = {int(q): k for q, k in T['cls'].items()}
ports = T['ports']; v = 112
dot = lambda a, b: bin(a & b).count('1') & 1
def rank(vs):
    bb = []
    for y in vs:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return bb
srcs = {i: tuple(s for s in range(v) if len(rank(r['plane'] + [ports[s]])) == 3) for i, r in enumerate(R)}
par = {}
def find(a):
    while par.get(a, a) != a: a = par[a]
    return a
for S in srcs.values():
    for s in S[1:]: par[find(s)] = find(S[0])
blk = {i: find(S[0]) for i, S in srcs.items()}
Q = {}
for i, S in srcs.items(): Q.setdefault(blk[i], set()).update(ports[s] for s in S)
Q = {b: rank(sorted(x)) for b, x in Q.items()}
planes_of = defaultdict(set)
for i, r in enumerate(R):
    for q, k, tv in r['visitors']:
        for (t, ph, role, p, f) in ev[q]:
            if t > tv and role == 'r' and cls[p] == 'y' and 11003 <= t < 12921: planes_of[p].add(i)
inperp = lambda F, t: all(not dot(x, ports[t - v]) for x in F)
c1 = Counter(); c2 = Counter(); c3 = Counter()
for t, s in planes_of.items():
    bs = {blk[i] for i in s}
    c1[len(bs)] += 1
    for b in bs: c2[inperp(Q[b], t)] += 1
for b in Q: c3[sum(1 for t in range(v, 2*v) if inperp(Q[b], t))] += 1
print('targets by number of distinct blocks among their B-consumed packets:', dict(c1))
print('(target, block) consumption pairs with Q_b inside u_t^perp:', dict(c2))
print('blocks by number of targets whose perp contains Q_b:', dict(c3))
# rank-5 relaxation: the smallest frames containing P that lie in all consumers' perps
c4 = Counter()
for i, r in enumerate(R):
    ts = [t for t in planes_of if i in planes_of[t]]
    if not ts: continue
    inter = [x for x in range(1, 512) if all(not dot(x, ports[t - v]) for t in ts)]
    c4[(len(ts), len(rank(inter)))] += 1
print('(consumers of P, dim of the intersection of their perps) :', dict(c4))
