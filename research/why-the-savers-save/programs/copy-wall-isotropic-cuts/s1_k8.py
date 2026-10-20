"""s1_k8.py -- (a) exact clique partitions of C_h (fewest orthogonal cliques), (b) the K = 8 survivors of clique counting
in model M1 (full stars + full line-cut totals, chain rule): need omega(P112 minus the classed ports) <= K0."""
from base import *
from collections import Counter
from itertools import combinations
TI = ti_subspaces(1)
iso = [b[0] for b in TI[1].values()]
Cs = {h: set(C_of([h])) for h in iso}
coord = [UI[o] for o in range(8)]
print('removed-octad labels o_k:', coord, ' isotropic?', [not dot9(o, o) for o in coord])
special = [h for h in iso if len(Cs[h]) == 56 and omega(Cs[h]) == 7]
print('the 8 classes with omega 7:', special)
for h in special:
    comp = set(range(N)) - Cs[h]
    rel = [k for k in range(8) if set(i for i in range(N) if not dot9(lab[i], coord[k])) == Cs[h]]
    print('  h %3d  weight %d  |C_h| %d  omega(complement) %d  C_h equals C_{o_k} for k in %s  h.o_k = %s' %
          (h, bin(h).count('1'), len(Cs[h]), omega(comp), rel, [dot9(h, o) for o in coord]))
# K0 = 7, one class: need omega(P112 \ C_h) <= 7
print('(K0=7, one class) omega of complement of C_h, by |C_h|:')
cnt = Counter((len(Cs[h]), omega(set(range(N)) - Cs[h])) for h in iso)
for k, c in sorted(cnt.items()): print('   |C_h| %d  omega(complement) %d : %d' % (k + (c,)))
# K0 = 6, two classes: need omega(P112 \ (C_h1 u C_h2)) <= 6
best = Counter()
for h1, h2 in combinations(iso, 2):
    rest = set(range(N)) - Cs[h1] - Cs[h2]
    best[omega(rest)] += 1
print('(K0=6, two classes) omega of the unclassed rest over all pairs:', sorted(best.items()))
