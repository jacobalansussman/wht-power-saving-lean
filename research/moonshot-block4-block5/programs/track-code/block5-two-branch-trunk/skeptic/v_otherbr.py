"""v_otherbr.py -- skeptic v-b5-twobranch (own test logic; F2 primitives from the claimant's sp2.py / fastP.py, copied to lib/).
Question: on the lam1 = 4 CLOSE branch N_c, can a one-array companion inside S_Z arrive through a THIRD branch j (frame
N_j != N_c)?  Such a companion is NOT multiplied by the same-frame twin factor G = mu_c + mu_c' u, so with a twin
(u outside S_Z) Lemma Z's proper-ideal argument (and Lemma Z' 'escape needs companions outside S_Z') does not apply.

Global coordinates (claimant's): t1 main at L0, s at M0 = Dbar L0, t2 at N1 = c2bar^-1 M0, branch k at
N_k = Theta_k(L0), Theta_k = c2bar^-1 Dbar c3kbar^-1 Dbar (Theta_k(M0) = N1).  Row of column x on a term with frame map Phi
sits at Phi(xbar^-1 M0); class direction Z = Theta_c^-1(A1).  x2 t1-injection: transport g with g(M0) = Q, all such g are
qbar P Dbar (close4_sp.py).  Its branch-j term sits at g(N_j) and has frame map g Theta_j.
Companion at N_c via j: g(N_j) = N_c.  The column x whose companion term hits the L-row at A1 has xbar^-1 M0 =
Theta_j^-1 g^-1 (A1); its branch-c array-1 main row is at Theta_c Theta_j^-1 g^-1(A1); S_Z <=> g^-1(A1) = W := Theta_j Theta_c^-1 (A1).
So: need b in P with bbar(qbar^-1 N_c) = Dbar N_j and bbar(qbar^-1 A1) = Dbar W   (g = qbar bbar^-1 Dbar).
Sanity: j = c gives W = A1, i.e. close4_sp.py's own-branch S_Z test.
usage: v_otherbr.py seed ncfg nb"""
import sys, os, random, time
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
from sp2 import mul, inv, image, dist, neighbors, rand_P, D, L0
from fastP import fast_find_P as fp

seed, ncfg, nb = [int(a) for a in sys.argv[1:4]]
rng = random.Random(seed); t0 = time.time()
M0 = image(D, L0)
tal = Counter(); twin = Counter(); done = 0; tries = 0
while done < ncfg and tries < 200 * ncfg:
    tries += 1
    c2 = rand_P(rng)
    N1 = image(inv(c2), M0)
    if dist(M0, N1) != 2:
        continue
    # admissible branches: random c3 letters, N_k = c2^-1 D c3^-1 M0 with d(L0,N_k)=2, d(M0,N_k)=3
    def theta(c3):
        return mul(inv(c2), mul(D, mul(inv(c3), D)))
    anchors = [X for X in neighbors(M0) if dist(X, L0) == 2 and dist(X, N1) == 3]
    P1, Q = rng.sample(anchors, 2)
    if dist(P1, Q) != 2:
        continue
    # close branch letter c3c
    c3c = None
    for _ in range(400):
        c = rand_P(rng); Th = theta(c); N = image(Th, L0)
        if dist(L0, N) == 2 and dist(M0, N) == 3 and dist(P1, N) == 4 and dist(Q, N) == 3:
            c3c = c; break
    if c3c is None:
        continue
    Thc = theta(c3c); Nc = image(Thc, L0)
    assert image(Thc, M0) == N1
    A1s = [A for A in neighbors(Nc) if dist(P1, A) == 5 and dist(Q, A) == 4]
    if not A1s:
        continue
    A1 = rng.choice(A1s)
    Z = image(inv(Thc), A1)
    assert dist(Z, L0) == 1
    a2 = fp(rng, M0, image(D, Q)); qb = mul(D, mul(a2, D)); assert image(qb, L0) == Q
    X1 = image(inv(qb), Nc); X2 = image(inv(qb), A1)
    done += 1
    # (a) twins: ubar in Stab(L0, M0); image of twin-family main row Theta_c(ubar Z)
    for _ in range(nb):
        u = fp(rng, M0, M0); assert image(u, L0) == L0 and image(u, M0) == M0
        Y = image(Thc, image(u, Z))
        twin[(Y == A1, dist(P1, Y))] += 1
    # (b) third branches j: random letters c3j giving an admissible frame N_j != N_c
    nj = 0
    for _ in range(300):
        if nj >= 6:
            break
        cj = rand_P(rng); Thj = theta(cj); Nj = image(Thj, L0)
        if not (dist(L0, Nj) == 2 and dist(M0, Nj) == 3):
            continue
        same = (Nj == Nc)
        W = image(mul(Thj, inv(Thc)), A1)
        Y1 = image(D, Nj); Y2 = image(D, W)
        key0 = ('same' if same else 'other', dist(P1, Nj), dist(Q, Nj), dist(M0, W))
        hit = 0; ok = 0
        for _ in range(nb):
            b = fp(rng, X1, Y1)
            if b is None:
                break
            ok += 1
            g = mul(qb, mul(inv(b), D))
            assert image(g, M0) == Q and image(g, Nj) == Nc
            if image(b, X2) == Y2:
                hit += 1
                assert image(inv(g), A1) == W
        nj += 1
        tal[key0 + ('tau found' if ok else 'no tau', 'S_Z hit' if hit else 'no S_Z hit')] += 1
        if hit:
            tal[('hits/samples',) + key0] += 0  # marker
            tal[('S_Z-rate', key0)] += hit
            tal[('samples', key0)] += ok
    # sanity: own branch j = c (close4 analogue): W = A1
    hit = 0
    for _ in range(nb):
        b = fp(rng, X1, image(D, Nc))
        if image(b, X2) == image(D, A1):
            hit += 1
    tal[('own-branch S_Z hits', nb)] += hit
    print("cfg %d done (%.0fs)" % (done, time.time() - t0), flush=True)
print("configs", done, "tries", tries)
print("twins (Y == A1 i.e. u in S_Z, d(P1, Theta_c(ubar Z))):", dict(sorted(twin.items())))
for k, v in sorted(tal.items(), key=str):
    print("  ", k, v)
