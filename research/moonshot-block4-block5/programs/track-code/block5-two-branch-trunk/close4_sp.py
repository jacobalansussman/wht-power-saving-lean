"""close4_sp.py (b5-twobranch): exact F2 test for the lam1 = 4 CLOSE branch N_c (d(P1,N_c)=4, d(Q,N_c)=3, A1 ~ N_c).
Lemma Z' (PROOF sec. 2b): within the class coset, family c escapes Lemma Z only through companions w (array-2 t3 terms in
B h_c q) OUTSIDE S_Z whose producing columns y_{c beta w^-1} have their array-1 image Theta_c(wbar Z) at layer <= 4.
Own-branch companions (x2 t1-injection, tau(M0) = Q, tau(N_c) = N_c): image frame Y = tau^-1(A1), a neighbour of N_c with
d(M0, Y) = d(Q, A1) = 4.  Census of the neighbours Y of N_c with d(M0, Y) = 4 by (d(P1, Y), Y == A1), and of the frames Y
actually reached by sampled tau.  usage: close4_sp.py nsamples ntau seed"""
import sys, random, time
from collections import Counter
from sp2 import *
from fastP import fast_find_P as fp
ns, ntau, seed = [int(a) for a in sys.argv[1:4]]
rng = random.Random(seed); t0 = time.time()
M0 = image(D, L0)
N1s = [N for N in sorted(neighbors(L0)) if dist(M0, N) == 2]
cand = Counter(); reach = Counter(); done = 0
while done < ns:
    N1 = rng.choice(N1s)
    brs = [N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3]
    anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3]
    P1, Q = rng.sample(anchors, 2)
    if dist(P1, Q) != 2:
        continue
    closes = [N for N in brs if dist(P1, N) == 4 and dist(Q, N) == 3]
    if not closes:
        continue
    Nc = rng.choice(closes)
    nb = sorted(neighbors(Nc))
    A1s = [A for A in nb if dist(P1, A) == 5 and dist(Q, A) == 4]
    if not A1s:
        continue
    A1 = rng.choice(A1s)
    for Y in nb:
        if dist(M0, Y) == 4:
            cand[(dist(P1, Y), Y == A1)] += 1
    # tau with tau(M0) = Q, tau(Nc) = Nc: tau = pibar beta^-1 Dbar; equivalently tau = g with g M0 = Q, g Nc = Nc.
    # sample: tau = s2 h s1 where s1: M0 -> Q (fixed choice), then h in Stab(Q) with h(s1 Nc) = Nc
    found = 0
    a2 = fp(rng, M0, image(D, Q))
    qb = mul(D, mul(a2, D))                     # qbar with qbar L0 = Q
    assert image(qb, L0) == Q
    X = image(inv(qb), Nc); Yt = image(D, Nc)
    for _ in range(ntau):
        b = fp(rng, X, Yt)                      # bbar X = Dbar Nc  ->  tau = qbar bbar^-1 Dbar fixes Nc, tau(M0) = Q
        if b is None:
            break
        g = mul(qb, mul(inv(b), D))
        assert image(g, M0) == Q and image(g, Nc) == Nc
        found += 1
        Y = image(inv(g), A1)
        reach[(dist(P1, Y), Y == A1, dist(M0, Y))] += 1
    done += 1
    print("sample %d: tau found %d of %d (%.0fs)" % (done, found, ntau, time.time() - t0), flush=True)
print("neighbours Y of N_c with d(M0,Y)=4, by (d(P1,Y), Y == A1):", dict(sorted(cand.items())))
print("frames Y = tau^-1(A1) reached by sampled tau, by (d(P1,Y), Y == A1, d(M0,Y)):", dict(sorted(reach.items())))
