"""split_sp.py (b5-twobranch): EXACT F_2 census of the SPLIT-TRUNK geometry for the engine 2 in 7 at block 5.

Shape: s_l = D sig_l x_l; t1 = D(a1 s1 + a2 s2 + inj); t2 = D c2 t1; t3 = D(c3a + mu c3b) t2; ends e1, e2.
Global picture (w3-b5 PROOF section 1): s at M0, t1 at L0, t2 at N1, the two t3 branches at Na = Ha L0, Nb = Hb L0
(Ha = c2bar^-1 Dbar c3abar^-1 Dbar = the global map Theta_a of branch a).  Array l's K sits at A_l (opposite anchor_l).
Joint progress is needed only up to N1 (both anchors at distance 3 from N1).  Array 1's class letters come from a branch
with d(P1, N) = 4 and A1 ~ N; array 2's from a branch with d(Q, N) = 4 and A2 ~ N.
Top-layer cleanliness of a class family of array 1 on branch a: Z = Ha^-1 A1 (a neighbour of L0); the other branch's
image Y = Hb Z must satisfy: d(P1, Y) <= 4 or Y = A1, and d(Q, Y) <= 4 or Y = A2 (else a layer-5 row with one producer).
usage: split_sp.py nwalks nletters seed"""
import sys, random, time
from collections import Counter
from sp2 import *
from fastP import fast_find_P


def main():
    nw, nlet, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(seed); t0 = time.time()
    M0 = image(D, L0)
    nbL0 = sorted(neighbors(L0))
    N1s = [N for N in nbL0 if dist(M0, N) == 2]
    prof = Counter(); conf = Counter(); top = Counter(); examples = []
    for w in range(nw):
        N1 = rng.choice(N1s)
        brs = [N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3]
        anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3]
        c2 = fast_find_P(rng, N1, M0)
        Hs = {}
        for N in brs:
            Hn = image(D, image(c2, N))
            lst = []
            for _ in range(nlet):
                c3 = fast_find_P(rng, Hn, M0)
                H = mul(inv(c2), mul(D, mul(inv(c3), D)))
                assert image(H, L0) == N
                lst.append(H)
            Hs[N] = lst
        nb = {N: sorted(neighbors(N)) for N in brs}
        for _ in range(40):
            P1, Q = rng.sample(anchors, 2)
            dPQ = dist(P1, Q)
            pr = {N: (dist(P1, N), dist(Q, N)) for N in brs}
            for N in brs:
                prof[(dPQ,) + pr[N]] += 1
            # branch pairs (Na for array 1, Nb for array 2), Na != Nb
            for Na in brs:
                if pr[Na][0] != 4:
                    continue
                A1s = [A for A in nb[Na] if dist(P1, A) == 5 and dist(Q, A) <= 4]
                if not A1s:
                    continue
                for Nb in brs:
                    if Nb == Na or pr[Nb][1] != 4:
                        continue
                    A2s = [A for A in nb[Nb] if dist(Q, A) == 5 and dist(P1, A) <= 4]
                    if not A2s:
                        continue
                    A1 = rng.choice(A1s); A2 = rng.choice(A2s)
                    lam = (dist(Q, A1), dist(P1, A2))
                    key = (dPQ, pr[Na], pr[Nb], lam)
                    conf[key] += 1
                    # top-layer cleanliness over sampled letter pairs
                    for i in range(nlet):
                        Ha = Hs[Na][i]; Hb = Hs[Nb][(i * 7 + 3) % nlet]
                        Za = image(inv(Ha), A1); Ya = image(Hb, Za)          # array-1 class family seen by branch b
                        Zb = image(inv(Hb), A2); Yb = image(Ha, Zb)          # array-2 class family seen by branch a
                        ok1 = (dist(P1, Ya) <= 4 or Ya == A1) and (dist(Q, Ya) <= 4 or Ya == A2)
                        ok2 = (dist(P1, Yb) <= 4 or Yb == A1) and (dist(Q, Yb) <= 4 or Yb == A2)
                        top[key + (ok1 and ok2, (dist(P1, Ya), dist(Q, Ya)), (dist(P1, Yb), dist(Q, Yb)))] += 1
                        if ok1 and ok2 and len(examples) < 5:
                            examples.append((w, key, (dist(P1, Ya), dist(Q, Ya)), (dist(P1, Yb), dist(Q, Yb))))
        print("walk %d done (%.0fs): %d branch frames, %d anchors" % (w, time.time() - t0, len(brs), len(anchors)), flush=True)
    print("branch profiles (d(P1,Q), d(P1,N), d(Q,N)) over sampled anchor pairs:")
    for k, v in sorted(prof.items()):
        print("   ", k, v)
    print("split configurations (d(P1,Q), prof(Na), prof(Nb), lam=(d(Q,A1), d(P1,A2))):")
    for k, v in sorted(conf.items()):
        print("   ", k, v)
    tot = Counter(); okc = Counter()
    for k, v in top.items():
        tot[k[:4]] += v
        if k[4]:
            okc[k[:4]] += v
    print("top-layer clean fraction by configuration type:")
    for k in sorted(tot):
        print("   ", k, "%d / %d" % (okc[k], tot[k]))
    dist_img = Counter()
    for k, v in top.items():
        dist_img[(k[5], k[6])] += v
    print("other-branch image layers ((d(P1,Ya), d(Q,Ya)), (d(P1,Yb), d(Q,Yb))):", dict(sorted(dist_img.items())))
    print("examples:", examples)
    print("time %.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
