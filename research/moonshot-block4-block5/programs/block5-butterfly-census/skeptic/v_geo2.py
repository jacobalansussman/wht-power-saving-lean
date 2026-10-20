"""v_geo2.py: check census claim 'ordered anchor pair admits a special-position merge iff P1 cap Q cap N2 = 0' on 60 N2s,
plus that the L3 fact's Y-set (nbr N2, d(M0,Y)=4) has 48 elements and the (d(P1,Y), d(Q,Y)) distribution."""
import sys
from collections import Counter
sys.path.insert(0, '.')
from v_geo import Lag, nbrs, dist, bc, m
L0 = Lag([1 << i for i in range(m)])
nL0 = nbrs(L0); M0 = nL0[7]; nM0 = nbrs(M0)
N2s = {}
for N1 in nL0:
    if dist(M0, N1) != 2: continue
    for N2 in nbrs(N1):
        if dist(L0, N2) == 2 and dist(M0, N2) == 3: N2s.setdefault(N2.mask, N2)
tab = Counter(); dist_tab = Counter()
for idx, N2 in enumerate(list(N2s.values())[::19]):
    nN2 = nbrs(N2); anchors = [P for P in nM0 if dist(P, N2) == 4]
    far = [Y for Y in nN2 if dist(M0, Y) == 4]
    for P in anchors:
        for Q in anchors:
            if P is Q: continue
            A1 = [Y for Y in nN2 if dist(P, Y) == 5 and dist(Q, Y) <= 4]
            A2 = [Y for Y in nN2 if dist(Q, Y) == 5 and dist(P, Y) <= 4]
            special = bool(A1) and bool(A2)
            triv = bc(P.mask & Q.mask & N2.mask) == 1
            tab[(special, triv)] += 1
            if special:
                for Y in far: dist_tab[(dist(P, Y), dist(Q, Y))] += 1
print("N2 sampled:", len(list(N2s.values())[::19]))
print("(special, P1capQcapN2==0) ->", dict(tab))
print("far-Y (d(P1,Y), d(Q,Y)) over special pairs:", dict(dist_tab))
