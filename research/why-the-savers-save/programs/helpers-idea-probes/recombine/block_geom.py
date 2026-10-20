#!/usr/bin/env python3
"""g-recombine: geometry of the Hamming blocks of the U8 word (planes = the 168 rank-3 frames the word visits, as in
u-helpers-hub/toy_visits.py), for the pair-factory (mutable sources) and butterfly / super-frame ideas.
Usage: python3 -I block_geom.py CANONICAL-WORD.json.gz"""
import sys, json, gzip, itertools
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
dot = lambda a, b: bin(a & b).count('1') & 1
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
planes = sorted({tuple(B[e[-1]]) for ph in ['A', 'B', 'tail', 'cleanup', 'cut', 'compensation', 'precondition']
                 for e in (c.get(ph) or []) if e[-1] is not None and len(B[e[-1]]) == 3
                 and sum(1 for s in range(v) if inspan(ports[s], list(B[e[-1]]))) == 4})
print('rank-3 gate frames holding exactly 4 labels:', len(planes))
srcs = {P: tuple(s for s in range(v) if inspan(ports[s], list(P))) for P in planes}
par = list(range(v))
def find(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
for P in planes:
    for a in srcs[P][1:]: par[find(a)] = find(srcs[P][0])
blocks = defaultdict(list)
for P in planes: blocks[find(srcs[P][0])].append(P)
print('blocks:', Counter(len(b) for b in blocks.values()))
cons = lambda S: [t for t in range(v) if all(not dot(x, ports[t]) for x in S)]
tot = Counter()
for bi, (k, Ps) in enumerate(sorted(blocks.items())):
    S = sorted({s for P in Ps for s in srcs[P]}); labs = [ports[s] for s in S]
    best = 0
    for d in {ports[a] ^ ports[b] for a, b in itertools.combinations(S, 2)}:
        pairs = [(a, b) for a, b in itertools.combinations(S, 2) if ports[a] ^ ports[b] == d]
        cov = sum(1 for P in Ps if any(set(p) | set(q) == set(srcs[P]) for p, q in itertools.combinations(pairs, 2)))
        best = max(best, cov)
    # can a pairing (perfect matching of the 8 sources, any differences) make every plane a union of two pairs?
    def matchings(L):
        if not L: yield []; return
        a = L[0]
        for i in range(1, len(L)):
            for m in matchings(L[1:i] + L[i+1:]): yield [(a, L[i])] + m
    bestm = max(sum(1 for P in Ps if any(set(p) | set(q) == set(srcs[P]) for p, q in itertools.combinations(m, 2))) for m in matchings(S))
    j = Counter()
    for P1, P2 in itertools.combinations(Ps, 2):
        U = list(P1) + list(P2); j[(len(basis(U)), len(cons(U)))] += 1
    if bi < 2:
        print(f'block {bi}: 8 sources span rank {len(basis(labs))}; planes {len(Ps)}; consumers per plane {Counter(len(cons(P)) for P in Ps)}')
        print(f'   best single difference class covers {best} planes as 2-pair sums; best perfect pairing covers {bestm}')
        print(f'   plane pairs (join rank, consumers of join): {dict(sorted(j.items()))}; consumers of block span {len(cons(labs))}')
    tot[(len(basis(labs)), best, bestm, len(cons(labs)))] += 1
print('all blocks (span rank, best d-class cover, best pairing cover, block-span consumers):', dict(tot))
