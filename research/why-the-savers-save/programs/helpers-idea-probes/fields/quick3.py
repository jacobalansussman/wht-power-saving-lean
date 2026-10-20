#!/usr/bin/env python3
"""g-fields quick screens: (a) register-allocation lens: live-range overlap of the 504 helpers (max clique of the interval
interference graph = registers needed if a register could be reused after its last event, frames ignored); live ranges
cut at the end of B (delivery) to see what cleanup costs; (b) CSE lens: pairs of packet planes read together by >= 2
targets in B; (c) sparse-factorisation lens: rank of each Hamming block's 8 labels and how many 4-label planes it has.
Usage: python3 -I quick3.py"""
import json, itertools
from collections import Counter, defaultdict
T = json.load(open('timeline.json')); R = json.load(open('consumers.json'))
ev = {int(q): E for q, E in T['ev'].items()}; cls = {int(q): k for q, k in T['cls'].items()}
ports = T['ports']; v = 112
H = [q for q in cls if q >= 2*v]
def maxover(iv):
    pts = sorted([(a, 1) for a, b in iv] + [(b + 0.5, -1) for a, b in iv]); m = c = 0
    for _, d in pts: c += d; m = max(m, c)
    return m
iv = [(ev[q][0][0], ev[q][-1][0]) for q in H if ev[q]]
print('(a) helpers', len(H), 'max simultaneously live (whole word):', maxover(iv))
ivB = [(a, min(b, 12920)) for a, b in iv if a <= 12920]
print('    if every helper were restored right after its last use before cleanup (cleanup ignored):', maxover([(ev[q][0][0], max([t for t, *_ in ev[q] if t < 13092] or [ev[q][0][0]])) for q in H]))
# by phase windows
for name, (a, b) in dict(A=(0, 4283), comp=(4683, 11003), B=(11003, 12921), cleanup=(13092, 19197)).items():
    act = sum(1 for q in H if any(a <= t < b for t, *_ in ev[q]))
    print('    helpers with an event in', name, act)
# (b) CSE
planes_of = defaultdict(set)
for i, r in enumerate(R):
    for q, k, tv in r['visitors']:
        for (t, ph, role, p, f) in ev[q]:
            if t > tv and role == 'r' and cls[p] == 'y' and 11003 <= t < 12921: planes_of[p].add(i)
pc = Counter()
for p, s in planes_of.items():
    for a, b in itertools.combinations(sorted(s), 2): pc[(a, b)] += 1
print('(b) plane pairs co-read by targets in B: by multiplicity', dict(Counter(pc.values())))
def rank(vs):
    bb = []
    for y in vs:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return len(bb)
shared = Counter()
for (a, b), m in pc.items():
    if m >= 2: shared[rank(R[a]['plane'] + R[b]['plane'])] += 1
print('    rank of P1+P2 for pairs co-read by >= 2 targets:', dict(shared))
# (c) blocks
lab = {}
for i, r in enumerate(R):
    S = [s for s in range(v) if rank(r['plane'] + [ports[s]]) == 3]
    lab[i] = tuple(S)
blocks = defaultdict(set)
for i, S in lab.items(): blocks[min(S)].add(i)
par = {}
def find(a):
    while par.get(a, a) != a: a = par[a]
    return a
for S in lab.values():
    for s in S[1:]: par[find(s)] = find(S[0])
bl = defaultdict(set)
for S in lab.values():
    for s in S: bl[find(s)].add(s)
rk = Counter(rank([ports[s] for s in B]) for B in bl.values())
print('(c) blocks', len(bl), 'sources per block', Counter(len(B) for B in bl.values()), 'label rank per block', dict(rk))
B0 = sorted(next(iter(bl.values())))
L = [ports[s] for s in B0]
n4 = sum(1 for c4 in itertools.combinations(range(8), 4) if rank([L[i] for i in c4]) == 3)
print('    block 0 labels', L, '4-subsets of rank 3 (candidate face planes):', n4,
      ' affine? sums of 4-subsets that are 0:', sum(1 for c4 in itertools.combinations(L, 4) if c4[0] ^ c4[1] ^ c4[2] ^ c4[3] == 0))
