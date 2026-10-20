import sys, time, collections
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-bit-side'))
import numpy as np
from e8geo import *
P = PUNCT_S2; live = [i for i in range(NL) if i not in P]
C = cliques(live, 8); print('cliques', len(C))
inc = np.zeros((len(C), NL), bool)
for i, c in enumerate(C): inc[i, list(c)] = True
t0 = time.time(); F = flats(3, verbose=True)
for r in (1, 2, 3):
    best = 0
    for m in F[r]:
        b = bits(m); met = int(inc[:, b].any(1).sum()); best = max(best, met)
    print('rank', r, 'flats', len(F[r]), 'most cliques met', best, '%.1f s' % (time.time() - t0), flush=True)
