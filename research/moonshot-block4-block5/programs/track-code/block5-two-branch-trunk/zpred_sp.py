"""zpred_sp.py (b5-twobranch): predecessor frames of the main rows of z_{t3,b} columns (b not in B_D').
Row b h_k pi (array l) has Pred = Theta_k(Z_b), Z_b = bbar^-1 M0 (a neighbour of L0, != M0), i.e. a neighbour Y of N_k
other than N1.  It has a y-producer only through an array-l term of layer >= layer(b h_k pi) - 1 at Y; the row is at layer
d(anchor_l, N_k), so for d(anchor_l, N_k) = 4 we need d(anchor_l, Y) = 3.  Census over all branch frames N_k (admissible
neighbours of N1) by profile (d(P1,N_k), d(Q,N_k)): number of neighbours Y != N1 with the needed distances for BOTH arrays."""
import sys, random, time
from collections import Counter
from sp2 import *
ns, seed = int(sys.argv[1]), int(sys.argv[2])
rng = random.Random(seed); t0 = time.time()
M0 = image(D, L0)
N1s = [N for N in sorted(neighbors(L0)) if dist(M0, N) == 2]
tab = Counter(); prof = Counter()
for s in range(ns):
    N1 = rng.choice(N1s)
    brs = [N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3]
    anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3]
    while True:
        P1, Q = rng.sample(anchors, 2)
        if dist(P1, Q) == 2:
            break
    for N in brs:
        pr = (dist(P1, N), dist(Q, N)); prof[pr] += 1
        ok = 0; tot = 0
        for Y in neighbors(N):
            if Y == N1:
                continue
            tot += 1
            need1 = (pr[0] < 4) or dist(P1, Y) == 3      # array-1 row at layer 4 needs a layer-3 array-1 term at Y
            need2 = (pr[1] < 4) or dist(Q, Y) == 3
            if need1 and need2:
                ok += 1
        tab[(pr, ok)] += 1
    print("sample %d (%.0fs)" % (s, time.time() - t0), flush=True)
print("branch profiles:", dict(sorted(prof.items())))
print("(profile, #neighbours Y != N1 giving both main rows a possible y-producer frame):")
for k, v in sorted(tab.items()):
    print("  ", k, v)
