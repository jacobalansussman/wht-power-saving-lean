# L-orbits on each (dist to L_Z, dist to L_X) class of Lagrangians: is the singular 7 the only L-orbit of size 7?
import itertools
from collections import Counter
from lib514 import *
lag = set()
for vs in itertools.combinations(range(1, NV), 3):
    S = span(vs)
    if len(S) == 8 and all(BTAB[a, b] == 0 for a in S for b in S): lag.add(S)
def dist(F1, F2): return 3 - {1: 0, 2: 1, 4: 2, 8: 3}[len(F1 & F2)]
cls = {}
for F in lag: cls.setdefault((dist(LZ, F), dist(LX, F)), []).append(F)
for c in sorted(cls):
    Fs = set(cls[c]); orbs = []
    while Fs:
        F = Fs.pop(); o = {frozenset(int(l[v]) for v in F) for l in LEVI}; Fs -= o; orbs.append(o)
    print(c, 'L-orbit sizes', sorted(len(o) for o in orbs),
          ' singular per orbit', [sum(all(q0(v) == 0 for v in F) for F in o) for o in sorted(orbs, key=len)])
