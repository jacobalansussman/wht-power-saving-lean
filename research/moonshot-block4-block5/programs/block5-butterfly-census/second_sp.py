"""second_sp.py (w3-b5): EXACT check of the second-order condition for lam = (3,3), over ALL frames (no sampling).
The first-order producer of L2 (x2 t1-injection, transport tau with tau M0 = Q, tau N2 = A1) brings rows b_j h p, b_j h q
whose predecessor frame is Y = tau^-1(N2) (NOTE section 3); Y is a neighbour of N2 with d(M0, Y) = d(Q, N2) = 4.  Those rows
have a second producer only if V has layer-3 terms of BOTH arrays at Y, i.e. d(P1, Y) = d(Q, Y) = 3.
Mirror (L1): x1 t1-injection with tau(M0) = P1, tau(N2) = A2, same Y-set condition.
Tabulates, over every walk sampled, every good anchor pair, every Y in nbr(N2) with d(M0, Y) = 4: (d(P1,Y), d(Q,Y))."""
import sys, random
from collections import Counter
from sp2 import *

rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
M0 = image(D, L0)
nw = int(sys.argv[1]) if len(sys.argv) > 1 else 3
tot = Counter(); totlow = Counter()
for w in range(nw):
    N1 = rng.choice([N for N in sorted(neighbors(L0)) if dist(M0, N) == 2])
    N2 = rng.choice([N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3])
    anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3 and dist(X, N2) == 4]
    nN2 = sorted(neighbors(N2))
    far = [Y for Y in nN2 if dist(M0, Y) == 4]
    dd = {(i, j): dist(X, Y) for i, X in enumerate(anchors) for j, Y in enumerate(nN2)}
    for i in range(len(anchors)):
        for k in range(len(anchors)):
            if i == k:
                continue
            A1s = [j for j in range(len(nN2)) if dd[(i, j)] == 5 and dd[(k, j)] == 3]
            A2s = [j for j in range(len(nN2)) if dd[(k, j)] == 5 and dd[(i, j)] == 3]
            if not (A1s and A2s):
                continue
            for j, Y in enumerate(nN2):
                if dist(M0, Y) == 4:
                    tot[(dd[(i, j)], dd[(k, j)])] += 1
                if dd[(i, j)] == 3 and dd[(k, j)] == 3:
                    totlow[dist(M0, Y)] += 1
    print("walk %d: |nbr(N2)| %d, |Y with d(M0,Y)=4| %d" % (w, len(nN2), len(far)), flush=True)
print("over lam=(3,3)-admissible ordered anchor pairs and all Y in nbr(N2) with d(M0,Y)=4: (d(P1,Y), d(Q,Y)) ->", dict(sorted(tot.items())))
print("low letters (d(P1,Y)=d(Q,Y)=3) by d(M0,Y):", dict(totlow))
