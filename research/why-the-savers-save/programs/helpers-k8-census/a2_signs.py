#!/usr/bin/env python3
"""g-search a2: per block plane, which targets need it and with which sign patterns (up to overall sign)."""
import sys, json
from fractions import Fraction as Fr
from collections import Counter, defaultdict
d = json.load(open('a0_demand.json')); ports = d['ports']; v = len(ports)
dem = {int(t): {int(s): Fr(x) for s, x in dd.items()} for t, dd in d['demand'].items()}
blocks = json.load(open('a1_blocks.json'))['blocks']; blk = {s: b for b, f in enumerate(blocks) for s in f}
dot = lambda a, b: bin(a & b).count('1') & 1
needs = defaultdict(list)  # (block, frozenset) -> list of (t, normalized sign tuple)
for t in range(v):
    per = defaultdict(dict)
    for s, x in dem[t].items(): per[blk[s]][s] = x
    for b, m in per.items():
        ss = sorted(m); x0 = m[ss[0]]
        pat = tuple(int(m[s]/x0) for s in ss)
        needs[(b, frozenset(ss))].append((t, pat, str(x0)))
sz = Counter(len(k[1]) for k in needs)
print('demanded block-subsets by size', sz)
for size in (4, 8):
    ks = [k for k in needs if len(k[1]) == size]
    print('size', size, 'subsets', len(ks), 'readers per subset', Counter(len(needs[k]) for k in ks))
    print('   distinct sign patterns per subset', Counter(len(set(p for _, p, _ in needs[k])) for k in ks))
    print('   readers per pattern', Counter(tuple(sorted(Counter(p for _, p, _ in needs[k]).values())) for k in ks))
# per block: planes used
pb = Counter(k[0] for k in needs if len(k[1]) == 4)
print('4-planes per block', Counter(pb.values()))
# are any 8-demands = union of two demanded 4-planes?
json.dump({f'{b}:{sorted(S)}': needs[(b, S)] for (b, S) in needs}, open('a2_needs.json', 'w'))
# label structure of a block: a + C
for b in range(2):
    L = [ports[s] for s in blocks[b]]
    print('block', b, 'sources', blocks[b], 'labels', [format(x, '09b') for x in L])
