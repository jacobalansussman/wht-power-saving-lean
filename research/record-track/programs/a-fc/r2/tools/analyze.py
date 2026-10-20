# fc-reopt round 2: structure of a #3945-line fill: per child, cap, window lengths, leftovers (old/new), patterns.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter, defaultdict
import numpy as np
import fillmodel as FM
M = FM.load(sys.argv[1]); fill = FM.read(sys.argv[2])['rows']
K = M['K']; ch = M['children']
used_o = defaultdict(int); used_n = defaultdict(int); pat = defaultdict(list)
for t, h, o, c, p, z in fill:
    b = (o, c, p); used_o[b] += int(M['W_old'][t, h]); used_n[b] += int(M['W_new'][t, h]); pat[b].append(int(M['W_new'][t, h]))
print('children (child: rank, family):', {c: (r, M['family'][c]) for c, r in sorted(ch.items())})
wl = defaultdict(Counter)
for (o, c), (cap, lo, hi) in M['windows'].items(): wl[c][hi - lo] += 1
for c in sorted(ch):
    lens = wl[c]; print('child %d rank %d window lengths: min %d max %d mean %.1f' % (c, ch[c], min(lens), max(lens), sum(k * v for k, v in lens.items()) / sum(lens.values())))
    left = Counter(); pats = Counter()
    for o in range(M['T']):
        for p in range(112):
            b = (o, c, p); left[(ch[c] - used_o.get(b, 0), ch[c] - used_n.get(b, 0))] += 1
            pats[tuple(sorted(pat.get(b, [])))] += 1
    print('   leftovers (old,new):', dict(sorted(left.items())))
    print('   top patterns (new widths):', pats.most_common(8))
