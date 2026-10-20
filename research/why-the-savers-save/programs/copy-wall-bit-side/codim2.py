"""codim2.py (d-bitspan): the root-orthogonal flats of E8: E7 hyperplanes l^perp (120) and the 6-flats
(l,m)^perp (E6 = A2^perp: 1120; D6 = A1A1^perp: 3780).  For each unit: live lines inside, cliques missed."""
import sys, collections, pickle; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import numpy as np
from e8geo import *
O = (G == 0)
def perp(lines):
    return [T for T in range(NL) if all(O[T, l] for l in lines)]
F7 = {}
for l in range(NL):
    F7[frozenset(perp([l]))] = ('E7', (l,))
F6 = {}
for l in range(NL):
    for m in range(l + 1, NL):
        f = frozenset(perp([l, m]))
        if f not in F6:
            F6[f] = ('E6' if G[l, m] != 0 else 'D6', (l, m))
print('E7 hyperplanes', len(F7), 'codim-2 6-flats', len(F6), collections.Counter(v[0] for v in F6.values()))
out = {}
for name, P in (('P111-s2', PUNCT_S2), ('P111-word', PUNCT_WORD), ('P120', [])):
    live = [i for i in range(NL) if i not in P]
    C = cliques(live, 8)
    Cs = [set(c) for c in C]
    res7 = {f: sum(1 for c in Cs if not (c & f)) for f in F7}
    res6 = {f: sum(1 for c in Cs if not (c & f)) for f in F6}
    print('==', name, 'cliques', len(C))
    print('  E7 hyperplanes: misses', sorted(collections.Counter(res7.values()).items()))
    for ty in ('E6', 'D6'):
        vals = [res6[f] for f in F6 if F6[f][0] == ty]
        print('  %s 6-flats: misses distribution (miss: count) %s' % (ty, sorted(collections.Counter(vals).items())[:12]))
        print('      min %d max %d mean %.1f  (fraction of cliques missed, min %.3f)' % (min(vals), max(vals), np.mean(vals), min(vals) / len(C)))
    # 6-flats defined by punctured lines
    if P:
        pp = [(a, b) for a in P for b in P if a < b]
        print('  6-flats (y_a, y_b)^perp from punctured pairs:', sorted((F6[frozenset(perp([a, b]))][0], int(G[a, b] != 0), res6[frozenset(perp([a, b]))]) for a, b in pp))
    out[name] = dict(res7={tuple(sorted(f)): v for f, v in res7.items()}, res6={tuple(sorted(f)): (F6[f], v) for f, v in res6.items()})
pickle.dump(out, open('out/codim2.pkl', 'wb'))
