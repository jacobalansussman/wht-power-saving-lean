"""km_fixed.py (t-km): flats of the P112 label matroid fixed by a prescribed automorphism g (port permutation from
gens.npz), by closure of unions of g-cycles.  usage: python -I km_fixed.py 7|3"""
import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from model import *
import numpy as np, itertools, json
from collections import Counter
p = int(sys.argv[1]); z = np.load('gens.npz'); perm = z['perm%d' % p]
cyc = []; seen = set()
for i in range(NP):
    if i in seen: continue
    c = []; j = i
    while j not in seen: seen.add(j); c.append(j); j = perm[j]
    cyc.append(c)
print('cycles', len(cyc), Counter(len(c) for c in cyc))
def closure(ports):
    b = f2basis([LAB[i] for i in ports])
    sp = {0}
    for x in b: sp |= {y ^ x for y in sp}
    return frozenset(i for i in range(NP) if LAB[i] in sp), len(b)
flats = {}
for m in range(1, 1 << len(cyc)):
    ports = [i for k, c in enumerate(cyc) if m >> k & 1 for i in c]
    cl, d = closure(ports)
    if len(cl) == len(ports): flats[cl] = d
print('g-fixed flats by dim:', sorted(Counter(flats.values()).items()))
json.dump([[sorted(f), d] for f, d in flats.items()], open('fixed_flats_g%d.json' % p, 'w'))
