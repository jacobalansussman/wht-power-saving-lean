#!/usr/bin/env python3
"""g-small: (D) packet frames at other ranks in P112 -- lines holding 3 port labels (u, v, u+v) and how they nest in
rank-3/rank-4 frames, i.e. whether one monotone chain could visit several packet frames of different ranks;
(F) 'jump-over' registers in the U8 word: registers whose chain holds a nonzero subspace of a packet plane P and
then a frame strictly above P without ever sitting at P (candidates to carry a proper-subspace residue so that a
second packet copy needs no third visitor at P).
Usage: python3 -I u8_nesting.py CANONICAL-WORD.json.gz"""
import sys, json, gzip, itertools
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
def inspan(x, bb):
    for z in bb: x = min(x, x ^ z)
    return x == 0
dot = lambda a, b: bin(a & b).count('1') & 1
P = set(ports)
print('label weights', Counter(bin(u).count('1') for u in ports), ' self-orthogonal labels', sum(1 for u in ports if not dot(u, u)))
deg = Counter(sum(1 for w in ports if w != u and not dot(u, w)) for u in ports); print('Gram-orthogonal readers per port', deg)
tri = {tuple(sorted((a, b, a ^ b))) for a, b in itertools.combinations(ports, 2) if (a ^ b) in P}
print('lines with 3 port labels:', len(tri))
orth_tri = [t for t in tri if not dot(t[0], t[1])]
print('  of which pairwise orthogonal (a reader of one can read all three):', len(orth_tri))
# targets t with the line inside u_t^perp
print('  targets whose perp contains such a line (histogram over lines):', Counter(sum(1 for u in ports if all(not dot(u, x) for x in t)) for t in tri).most_common(6))
# replay frames
B = [basis(f) for f in F]; dims = [len(b) for b in B]
def sub(a, b): return all(inspan(x, B[b]) for x in B[a])
cur = list(c['start']); chain = defaultdict(list)
for r in range(n): chain[r].append(cur[r])
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not sub(cur[r], f): return
    cur[r] = f; chain[r].append(f)
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        climb(e[1], e[-1])
        if e[0] == 'add': climb(e[2], e[-1])
for r in range(n): climb(r, c['final'][r])
planes = sorted({tuple(B[f]) for r in range(n) for f in chain[r] if dims[f] == 3})
pidx = {}
for r in range(n):
    for f in chain[r]:
        if dims[f] == 3: pidx.setdefault(tuple(B[f]), set()).add(r)
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
jump = Counter(); jump_planes = Counter(); jump_by_gap = Counter()
for Pl in planes:
    cnt = 0
    for r in range(n):
        if r in pidx[Pl]: continue
        fr = chain[r]
        below = [f for f in fr if dims[f] >= 1 and all(inspan(x, list(Pl)) for x in B[f])]
        above = [f for f in fr if dims[f] > 3 and all(inspan(x, B[f]) for x in Pl)]
        if below and above:
            jump[cls(r)] += 1; cnt += 1
            jump_by_gap[(dims[below[-1]], dims[above[0]])] += 1
    jump_planes[cnt] += 1
print('jump-over registers (subspace of P, then above P, never at P), by class, summed over planes:', dict(jump))
print('(last rank below, first rank above):', sorted(jump_by_gap.items()))
print('planes by number of jump-over registers:', sorted(jump_planes.items())[:20])
