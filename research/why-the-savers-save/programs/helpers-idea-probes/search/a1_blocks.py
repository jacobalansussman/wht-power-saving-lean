#!/usr/bin/env python3
"""g-search a1: block geometry of the P112 labels and the shape of every target read (by blocks / planes)."""
import sys, json
from collections import Counter, defaultdict
d = json.load(open(sys.argv[1])); ports = d['ports']; v = len(ports); dem = {int(k): x for k, x in d['demand'].items()}
dot = lambda a, b: bin(a & b).count('1') & 1
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
def rank(xs): return len(basis(xs))
# blocks: sources whose labels lie in a common coset a + C, rank(C)=3 ... find by: s~s' if they appear together in a 4-plane read?
# simpler: for each pair, differences; group by maximal affine 3-flats containing 8 labels
labels = set(ports); idx = {p: i for i, p in enumerate(ports)}
flats = set()
for i in range(v):
    for j in range(i+1, v):
        for k in range(j+1, v):
            a, b, cc = ports[i], ports[j], ports[k]
            C = basis([a ^ b, a ^ cc])
            if len(C) < 2: continue
            # extend to 3-flats
            for l in range(v):
                dd = ports[l] ^ a
                C3 = basis(C + [dd])
                if len(C3) != 3: continue
                # all 8 points of a + span(C3)
                sp = {0}
                for z in C3: sp |= {x ^ z for x in sp}
                pts = {a ^ x for x in sp}
                if pts <= labels: flats.add(frozenset(idx[p] for p in pts))
print('affine 3-flats fully inside the label set:', len(flats))
cov = Counter(s for f in flats for s in f)
print('flats per source', Counter(cov.values()))
# a partition into 14?
fl = sorted(flats, key=sorted)
# greedy exact cover by DFS
best = None
def dfs(chosen, used):
    global best
    if best: return
    if len(used) == v: best = list(chosen); return
    s = min(set(range(v)) - used)
    for f in fl:
        if s in f and not (f & used): dfs(chosen + [f], used | f)
dfs([], frozenset())
print('partition into flats found:', best is not None, len(best) if best else None)
blk = {s: bi for bi, f in enumerate(best) for s in f}
json.dump(dict(blocks=[sorted(f) for f in best]), open('a1_blocks.json', 'w'))
# demand by block: per target, per block, size
shape = Counter()
for t in range(v):
    per = Counter(blk[int(s)] for s in dem[t])
    shape[tuple(sorted(per.values()))] += 1
print('per-target demand block sizes:', shape)
own = Counter(len([s for s in dem[t] if blk[int(s)] == blk[t]]) for t in range(v))
print('demand inside the own block:', own)
# read shapes
rs = Counter()
for ph, t, r, c_, srcs, k in d['reads']:
    per = Counter(blk[s] for s in srcs)
    rs[(c_, tuple(sorted(per.values())))] += 1
for k_, x in sorted(rs.items(), key=lambda z: -z[1]): print(' read', k_, x)
