"""census_sp.py (w3-b5): EXACT F_2 census of special-position merges for the 2-in-7 butterfly at block 5.

Global picture (derivation in NOTE.md section 1): for array l with merge element pi_l (p or q), a term g has global frame
G(g) = pibar_l gbar^-1 L0, and its layer is d(anchor_l, G(g)) with anchor_l = pibar_l L0.  The shared moves sit at global
frames  s: M0 = Dbar L0,  t1: L0,  t2: N1,  t3: N2,  end letter b: Y (a neighbour of N2).  Array 1's target K x1 sits at
A_1 = pbar F*, array 2's at A_2 = qbar F*.  The top layer (mw2-b5beyond PROOF 3c) needs:
  geodesic walks anchor -> M0 -> L0 -> N1 -> N2 for both anchors (distances 1,2,3,4) [joint progress],
  A_l a neighbour of N2 opposite anchor_l (the clean letters exist),
  clean splits: d(Q, A_1) <= 4 and d(P1, A_2) <= 4.
The leftovers of the clean split: x2 at global A_1 (layer d(Q, A_1)), x1 at global A_2 (layer d(P1, A_2)).
usage: census_sp.py nwalks seed"""
import sys, random, time
from collections import Counter
from sp2 import *


def walk_census(N1, N2, M0, nbM0):
    anchors = [X for X in nbM0 if dist(X, L0) == 2 and dist(X, N1) == 3 and dist(X, N2) == 4]
    nN2 = sorted(neighbors(N2))
    dt = {(i, j): dist(X, Y) for i, X in enumerate(anchors) for j, Y in enumerate(nN2)}
    T = {i: [j for j in range(len(nN2)) if dt[(i, j)] == 5] for i in range(len(anchors))}
    nconf = 0; npairs = 0; lam = Counter(); per_pair = Counter()
    for i in range(len(anchors)):
        for k in range(len(anchors)):
            if i == k:
                continue
            A1s = [j for j in T[i] if dt[(k, j)] <= 4]
            A2s = [j for j in T[k] if dt[(i, j)] <= 4]
            c = len(A1s) * len(A2s)
            per_pair[(len(T[i]), len(A1s), len(A2s))] += 1
            if c:
                npairs += 1; nconf += c
                for a in A1s:
                    for b in A2s:
                        lam[(dt[(k, a)], dt[(i, b)])] += 1
    return len(anchors), Counter(len(T[i]) for i in T), npairs, nconf, lam, per_pair


def main():
    nw, seed = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed); t0 = time.time()
    M0 = image(D, L0)
    nbM0 = sorted(neighbors(M0)); nbL0 = sorted(neighbors(L0))
    N1s = [N for N in nbL0 if dist(M0, N) == 2]
    print("M0 = Dbar L0; |nbr(L0)| = %d, N1 candidates (d(M0, N1) = 2): %d" % (len(nbL0), len(N1s)))
    seen = Counter()
    for w in range(nw):
        N1 = rng.choice(N1s)
        N2s = [N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3]
        N2 = rng.choice(N2s)
        na, Tc, npairs, nconf, lam, per_pair = walk_census(N1, N2, M0, nbM0)
        key = (len(N2s), na, tuple(sorted(Tc.items())), npairs, nconf, tuple(sorted(lam.items())))
        seen[key] += 1
        if w < 2:
            print("walk %d: |N2 choices| %d, anchors %d, targets per anchor %s, ordered anchor pairs with a special-position "
                  "merge %d of %d, configurations (P1, Q, A1, A2) %d, leftover layers (d(Q,A1), d(P1,A2)) %s" % (
                      w, len(N2s), na, dict(Tc), npairs, na * (na - 1), nconf, dict(sorted(lam.items()))))
            print("   per ordered pair (|T(P1)|, #clean A1, #clean A2): %s" % dict(sorted(per_pair.items())))
    print("distinct census results over %d random walks: %d  (walk-independent iff 1)" % (nw, len(seen)))
    for k, v in seen.items():
        print("  %d walks: %s" % (v, k))
    print("time %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
