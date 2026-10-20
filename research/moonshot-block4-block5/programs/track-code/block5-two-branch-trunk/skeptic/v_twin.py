# v_twin.py -- skeptic (own code). Same-frame "twin" branches on the close branch N_c (lam1 = 4).
# A twin h_j = u h_c (u in B_D' = B cap D B D^-1, i.e. ubar in Stab(L0, M0)) puts its K-letters' main-term image at
# Theta_c(ubar Z), Z = Theta_c^-1(A1); as u ranges over Stab(L0,M0) these frames range over g(A1), g in Stab(N_c, N1)
# (Theta_c maps L0 -> N_c, M0 -> N1).  Here: neighbours X of N_c with d(N1, X) = 2 (the candidate orbit; superset if
# Stab(N_c,N1) is not transitive on them), X != A1, tallied by d(P1, X).  d(P1,X) <= 4 = twin escape open at first order.
import sys, random
from collections import Counter
from v_census import nbrs, dist, span_mask
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
rnd = random.Random(seed)
L0 = span_mask([1 << (5 + i) for i in range(5)])
M0 = rnd.choice(nbrs(L0)); N1 = rnd.choice([N for N in nbrs(L0) if dist(N, M0) == 2])
anchors = [P for P in nbrs(M0) if dist(P, L0) == 2 and dist(P, N1) == 3]
nN1 = nbrs(N1)
tal = Counter(); tal2 = Counter(); ncfg = 0
for P1 in anchors[:6]:
    branches = [N for N in nN1 if dist(P1, N) == 4]
    targets = set()
    for N in branches:
        for A in nbrs(N):
            if dist(P1, A) == 5: targets.add(A)
    targets = sorted(targets)[:40]
    for Qa in anchors:
        if Qa == P1 or dist(P1, Qa) != 2: continue
        for A in targets:
            if dist(Qa, A) != 4: continue
            C = [N for N in nN1 if dist(N, A) == 1]
            Nc = [N for N in C if dist(Qa, N) == 3]; assert len(Nc) == 1; Nc = Nc[0]
            Xs = [X for X in nbrs(Nc) if dist(N1, X) == 2 and X != A]
            assert len(Xs) == 59, len(Xs)
            ncfg += 1
            tal[tuple(sorted(Counter(dist(P1, X) for X in Xs).items()))] += 1
            tal2[tuple(sorted(Counter((dist(P1, X), dist(Qa, X)) for X in Xs).items()))] += 1
print(f"walk seed {seed}: {ncfg} lam1=4 configurations (P1, Q, A1) checked")
for k, v in sorted(tal.items()): print("  d(P1, X) histogram over the 59 X != A1:", k, v)
for k, v in sorted(tal2.items()): print("  (d(P1,X), d(Q,X)) histogram over the 59 X != A1:", k, v)
