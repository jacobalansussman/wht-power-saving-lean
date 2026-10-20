#!/usr/bin/env python3
"""d-pair-market: the Hamming blocks the U8 word ACTUALLY uses (sources grouped by the 168 rank-3 frames it visits,
as in u-helpers-hub/toy_visits.py), compared with g-search's a1 partition (a DFS exact cover of the labels by affine
3-flats), and the demand shape per block under both partitions.
Usage: python3 -I u8_blocks.py WORD  -> u8_blocks.json (same format as a1_blocks.json)"""
import sys, json, gzip, os
from collections import Counter, defaultdict
here = os.path.dirname(os.path.abspath(sys.argv[0]))
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; ports = c['ports']
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
def inspan(x, bb):
    for z in bb: x = min(x, x ^ z)
    return x == 0
B = [basis(f) for f in F]
planes = [tuple(b) for b in B if len(b) == 3]
srcs = {P: [s for s in range(v) if inspan(ports[s], list(P))] for P in planes}
assert all(len(x) == 4 for x in srcs.values())
parent = list(range(v))
def find(a):
    while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
    return a
for P in planes:
    for a in srcs[P][1:]: parent[find(a)] = find(srcs[P][0])
bl = defaultdict(list)
for s in range(v): bl[find(s)].append(s)
U8 = sorted(sorted(x) for x in bl.values())
A1 = json.load(open(here + '/a1_blocks.json'))['blocks']
print('U8 blocks:', len(U8), 'sizes', Counter(len(x) for x in U8))
print('same partition as g-search a1:', sorted(map(tuple, U8)) == sorted(map(tuple, A1)))
blkU = {s: i for i, x in enumerate(U8) for s in x}; blkA = {s: i for i, x in enumerate(A1) for s in x}
print('a1 blocks meeting each U8 block:', Counter(len({blkA[s] for s in x}) for x in U8))
# is every U8 block an affine 3-flat (8 labels, rank 4 span, XOR-closed in triples)?
def aff3(x):
    L = [ports[s] for s in x]; a = L[0]; C = basis([y ^ a for y in L])
    return len(C) == 3
print('U8 blocks that are affine 3-flats:', sum(aff3(x) for x in U8))
d = json.load(open(here + '/a0_demand.json')); dem = {int(t): {int(s) for s in dd} for t, dd in d['demand'].items()}
for name, blk in (('a1', blkA), ('U8', blkU)):
    sh = Counter()
    for t in range(v):
        per = Counter(blk[s] for s in dem[t]); sh[tuple(sorted(per.values()))] += 1
    print(name, 'per-target demand sizes by block:', sh.most_common(4))
# planes of U8 as 4-subsets: are they demanded 4-subsets under the U8 partition?
req = Counter()
for t in range(v):
    per = defaultdict(set)
    for s in dem[t]: per[blkU[s]].add(s)
    for b, S in per.items(): req[frozenset(S)] += 1
pl4 = {frozenset(x) for x in srcs.values()}
print('U8 planes that are demanded 4-subsets (U8 partition):', sum(1 for P in pl4 if P in req), 'of', len(pl4),
      ' demanded 4-subsets not among U8 planes:', sum(1 for S in req if len(S) == 4 and S not in pl4),
      ' demanded 8-subsets:', sum(1 for S in req if len(S) == 8))
json.dump(dict(blocks=U8), open(here + '/u8_blocks.json', 'w'))
