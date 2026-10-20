#!/usr/bin/env python3
"""geom.py (fp-capacity): exact geometry of what one paid move can serve, on the frame graph Gamma_m = DSp(2m,2).

(1) Lemma S (image-frame picture). A label with image frame F (level k = d(F, L0)) and a paid move with neighbour
    N of L0: the level becomes d(F, N).  Claim: over the 2(2^m - 1) neighbours N, exactly 2(2^m - 2^k) ascend,
    2^k - 1 are lateral and 2^k - 1 descend; ascent  <=>  F cap L0 not inside the hyperplane N cap L0.
(2) Joint ascent. A bundle of image frames F_1..F_j (levels k_i < m) is ascended by SOME move whenever
    sum_i (2^k_i - 1) < 2^m - 1 (union bound); three frames at level m - 1 can block every move.  Checked: the
    smallest blocking set of level-(m-1) frames.
(3) Progress profile (position picture, raw gauge, f = d(L0,.) - d(L1,.)): for each position X, the multiset of
    Delta f over its neighbours.  Only Delta f = +2 ("progress") and +1 carry weight in Theorem A.
(4) Rigid 2-bundles. Two live classes in one register keep a fixed relative isometry; with (A, A') the image of
    (L0, L1) under it, the bundle's moves contribute at most  rho = (max g - min g)/2,  g = (f + f_A)/2,
    f_A = d(A,.) - d(A',.), to Theorem A's bound; lam = the longest path that progresses for both (A,A') and
    (L0,L1).  Distribution over all ordered opposite pairs (m = 2, 3: exact; m = 4: sampled).
usage: python3 geom.py [m_max_exact] [samples_m4]"""
import sys, os, json, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "ledger"))
from cutconst import setup
import numpy as np
from collections import Counter

def run(m, sample=None, rng=None):
    labs, L0bits, D = setup(m)
    V = len(labs); i0 = labs.index(L0bits)
    opp0 = np.flatnonzero(D[i0] == m)
    nb0 = np.flatnonzero(D[i0] == 1)
    out = dict(m=m, frames=V, neighbours=len(nb0))
    # (1) Lemma S
    okS = True; prof = Counter()
    for F in range(V):
        k = int(D[i0, F]); d = D[F, nb0] - k
        a, z, dd = int((d == 1).sum()), int((d == 0).sum()), int((d == -1).sum())
        prof[(k, a, z, dd)] += 1
        if (a, z, dd) != (2 * (2 ** m - 2 ** k), 2 ** k - 1, 2 ** k - 1): okS = False
    out["lemmaS_ok"] = okS
    out["lemmaS_profile"] = {"level %d" % k: dict(frames=c, ascend=a, lateral=z, descend=dd) for (k, a, z, dd), c in sorted(prof.items())}
    # (2) smallest blocking set among level-(m-1) frames: every neighbour N fails to ascend at least one frame
    lev = np.flatnonzero(D[i0] == m - 1)
    asc = (D[np.ix_(lev, nb0)] == m)        # frame at level m-1 ascends under N iff d(F,N) = m
    best = None
    for j in ((1, 2, 3, 4) if m <= 3 else ()):
        for comb in itertools.combinations(range(len(lev)), j):
            if not asc[list(comb)].all(axis=0).any(): best = j; break
        if best: break
    out["min_blocking_set_level_m_minus_1"] = best
    # (3) progress profile over positions, raw gauge (L1 = an opposite of L0; all opposites are equivalent)
    i1 = int(opp0[0])
    f = D[i0] - D[i1]
    pp = Counter()
    for X in range(V):
        nb = np.flatnonzero(D[X] == 1)
        df = f[nb] - f[X]
        key = (int(D[i0, X]), int(D[i1, X])) + tuple(int((df == v).sum()) for v in (2, 1, 0, -1, -2))
        pp[key] += 1
    out["progress_profile"] = ["(d(L0,X), d(L1,X)) = (%d,%d): %d positions; neighbours with df=+2,+1,0,-1,-2: %s"
                               % (k[0], k[1], c, k[2:]) for k, c in sorted(pp.items())]
    # (4) rigid 2-bundles
    pairs = [(A, Ap) for A in range(V) for Ap in np.flatnonzero(D[A] == m)]
    if sample is not None and len(pairs) > sample:
        pairs = [pairs[i] for i in rng.sample(range(len(pairs)), sample)]
    adj = (D == 1)
    order = np.argsort(f)
    rhoC = Counter(); lamC = Counter(); fullpairs = 0
    for (A, Ap) in pairs:
        fA = D[A] - D[Ap]
        g = (f + fA) / 2.0
        rho = (g.max() - g.min()) / 2.0
        rhoC[rho] += 1
        # longest co-progress path: arcs X->Y with f(Y)=f(X)+2 and fA(Y)=fA(X)+2
        L = np.zeros(V, dtype=np.int64)
        # dynamic programming by increasing f
        best = 0
        for val in sorted(set(f.tolist())):
            Xs = np.flatnonzero(f == val)
            for X in Xs:
                nb = np.flatnonzero(adj[X] & (f == val + 2) & (fA == fA[X] + 2))
                for Y in nb:
                    if L[X] + 1 > L[Y]: L[Y] = L[X] + 1
        lam = int(L.max())
        lamC[lam] += 1
        if rho == m: fullpairs += 1
    n = len(pairs)
    out["bundle_pairs_examined"] = n
    out["bundle_rho_distribution"] = {str(k): v for k, v in sorted(rhoC.items())}
    out["bundle_lambda_distribution"] = {str(k): v for k, v in sorted(lamC.items())}
    out["bundle_mean_rho"] = sum(k * v for k, v in rhoC.items()) / n
    out["bundle_mean_lambda"] = sum(k * v for k, v in lamC.items()) / n
    out["bundle_share_full_rho_m"] = fullpairs / n
    return out

def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    s4 = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    rng = random.Random(1)
    res = []
    for m in range(2, mmax + 1):
        r = run(m); res.append(r); print(json.dumps(r, indent=1), flush=True)
    if s4:
        r = run(4, sample=s4, rng=rng); r["sampled"] = True; res.append(r); print(json.dumps(r, indent=1), flush=True)
    with open(os.path.join(HERE, "geom.json"), "w") as fh: json.dump(res, fh, indent=1)

if __name__ == "__main__":
    main()
