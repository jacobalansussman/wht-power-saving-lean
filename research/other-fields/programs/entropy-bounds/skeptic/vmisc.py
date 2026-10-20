# vmisc.py -- v-nc-gain: (i) the known gap layouts are not two-class stretchable; (ii) stretchable shares with more samples
import sys, random, collections
sys.path.insert(0, '.')
import numpy as np, vlemmas
from vcut_tests import Mapped
g2, g3 = Mapped(2), Mapped(3)
for g, P in [(g2, [(0,7),(3,1),(4,10),(0,10),(10,7)]), (g3, [(0,95),(0,80),(80,95),(84,60)]),
             (g3, [(0,95),(0,80),(80,95),(12,80),(84,60)]), (g3, [(0,128),(0,124),(124,128),(8,124),(76,44)]),
             (g3, [(4,123),(25,33),(78,48)]), (g3, [(34,15),(65,28),(126,113)]), (g3, [(114,97),(17,41),(16,29)])]:
    print(f"m={g.m} {P}: opposite {all(g.D[s,t]==g.m for s,t in P)}; two-class stretchable: {vlemmas.two_class(g, P) is not None}"
          + (f"; e(D) = {vlemmas.edef(g, P)[0]}" if len(P) == 3 else ''))
rng = random.Random(99)
for g, W, N in [(g3, 5, 2000), (g3, 4, 1000)]:
    hit = sum(1 for _ in range(N) if vlemmas.two_class(g, vlemmas.rand_layout(g, rng, W)) is not None)
    print(f"m={g.m} W={W}: {hit}/{N} = {100*hit/N:.1f}% two-class stretchable (distinct frames)")
