#!/usr/bin/env python3
"""g-break: frame annotation per phase; block spans and which target perps contain them."""
import sys, json, gzip
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
B = [basis(f) for f in F]; dims = [len(b) for b in B]
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
for ph in ['A', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    E = c.get(ph) or []
    fr = Counter('None' if e[-1] is None else dims[e[-1]] for e in E)
    kinds = Counter((e[0], cls(e[1]), cls(e[2]) if e[0]=='add' else '-') for e in E)
    print(ph, len(E), 'frame dims', dict(sorted(fr.items(), key=str)), 'kinds', dict(kinds))
print('frame ids >= 999 used in tail:', Counter(e[-1] for e in c['tail'])[999], len(F))
# blocks: 8 sources per block with rank-4 span
print('ports == target ports? targets use ports[t] (assumed)')
labs = ports
# components via planes containing 4 labels: recompute blocks via rank-4 span grouping
spans = defaultdict(list)
for s in range(v):
    pass
# group sources by the toy's blocks: use planes from census rule is slow; instead group by rank-4 span of 8 labels using the
# known block from explore.py: sources whose label lies in S
seen = set(); blocks = []
import itertools
for s in range(v):
    if s in seen: continue
    # find the 7 others: labels l such that some rank-3 subspace contains l, ports[s] and two more labels with no 3 on a line
    # simpler: block = sources with label in the rank-4 span generated greedily from planes in the word is complex; use pairwise
    pass
