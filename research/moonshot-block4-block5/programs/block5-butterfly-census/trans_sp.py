"""trans_sp.py (w3-b5): is Stab(L0, M0) transitive on the shared walks (N1, N2)?  (Then every census / second-order count
done for one walk holds for all 60 x 56 = 3360 walks: walk-independence becomes exact.)  BFS orbit of one walk under 40
random generators of Stab(L0) cap Stab(M0)."""
import random
from sp2 import *
rng = random.Random(11)
M0 = image(D, L0)
gens = []
while len(gens) < 40:
    g = rand_P(rng)
    if image(g, M0) == M0:
        gens.append(g)
N1s = [N for N in sorted(neighbors(L0)) if dist(M0, N) == 2]
walks = set()
for N1 in N1s:
    for N2 in neighbors(N1):
        if dist(L0, N2) == 2 and dist(M0, N2) == 3:
            walks.add((N1, N2))
start = next(iter(walks))
orb = {start}; fr = [start]
while fr:
    w = fr.pop()
    for g in gens:
        w2 = (image(g, w[0]), image(g, w[1]))
        if w2 not in orb:
            orb.add(w2); fr.append(w2)
print("walks M0-L0-N1-N2 (geodesic):", len(walks), " orbit of one walk under Stab(L0,M0):", len(orb), " transitive:", orb == walks)
