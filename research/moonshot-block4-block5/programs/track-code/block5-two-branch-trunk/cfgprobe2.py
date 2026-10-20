"""cfgprobe2.py: for random (walk, anchors), every A1 opposite P1 adjacent to an admissible branch (nbr of N1, d(L0)=2,
d(M0)=3): the profiles (d(P1,N), d(Q,N)) of ALL admissible branches adjacent to it, and lam1 = d(Q,A1)."""
import sys, random, time
from collections import Counter
from sp2 import *
seed, ntry = int(sys.argv[1]), int(sys.argv[2])
rng = random.Random(seed); t0 = time.time()
M0 = image(D, L0)
N1s = [N for N in sorted(neighbors(L0)) if dist(M0, N) == 2]
pat = Counter()
for t in range(ntry):
    N1 = rng.choice(N1s)
    brs = [N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3]
    anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3]
    while True:
        P1, Q = rng.sample(anchors, 2)
        if dist(P1, Q) == 2:
            break
    byA = {}
    for N in brs:
        if dist(P1, N) == 4:
            for A in neighbors(N):
                if dist(P1, A) == 5:
                    byA.setdefault(A, []).append(N)
    for A, Ns in byA.items():
        pat[(dist(Q, A), tuple(sorted(dist(Q, N) for N in Ns)))] += 1
    print("try %d: %d branches, %d targets (%.0fs)" % (t, len(brs), len(byA), time.time() - t0), flush=True)
for k, v in sorted(pat.items()):
    print("  lam1 = d(Q,A1) = %d, d(Q, N) of the branches adjacent to A1: %s : %d" % (k[0], k[1], v))
