"""even_test.py (ex2-families, session 2, 14 Oct 2026 23:05 EDT): the E-pattern tests (PROOF.md Theorem 10).

E = I + Gamma - K (M x M), K = Btilde Atilde (rank <= W), Gamma strictly lower triangular in time order.
Facts (PROOF.md section 10):
  (a) E vanishes at every position (i, j) with i, j in a common tight support (diagonal included);
  (b) even Q:  E[D_Q, D_Q] is NONSINGULAR;
  (c) odd Q:   rank E[D_Q, D_Q] = |D_Q| - W.
Pattern consequences tested here (order-free, coefficient-free):
  EVEN TEST: for every even Q with D_Q nonempty, the bipartite graph of positions of D_Q x D_Q not forced to 0 by (a)
             has a perfect matching (else every matrix with that pattern is singular, contradicting (b)).
  ODD TEST:  for every odd Q, that graph has a matching of size >= |D_Q| - W (term rank >= rank).
usage: even_test.py   (runs the calibration list and the block-4 survivors)
"""
import sys
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from refl import quad_forms
from clique_test import clique_test

_cache = {}


def forms(m):
    if m not in _cache:
        _cache[m] = quad_forms(m)
    return _cache[m]


def matching_size(allowed):
    if allowed.shape[0] == 0:
        return 0
    g = csr_matrix(allowed.astype(np.int8))
    mt = maximum_bipartite_matching(g, perm_type="column")
    return int((mt >= 0).sum())


def epattern_test(m, S, W, verbose=False, label=""):
    vals, tr = forms(m)
    sub = vals[:, S]
    wt = sub.sum(axis=1)
    odd = tr == 1
    M = len(S)
    dmin = int(wt[odd].min())
    if dmin < W:
        return ("deletion", dmin)
    tight = np.nonzero(odd & (wt == W))[0]
    zero = np.zeros((M, M), dtype=bool)
    for f in tight:
        D = np.nonzero(sub[f])[0]
        zero[np.ix_(D, D)] = True
    allowed = ~zero
    nE = nO = 0
    firstE = firstO = None
    for f in range(vals.shape[0]):
        D = np.nonzero(sub[f])[0]
        if len(D) == 0:
            continue
        ms = matching_size(allowed[np.ix_(D, D)])
        if tr[f] == 0:
            if ms < len(D):
                nE += 1
                if firstE is None:
                    firstE = (f, [S[i] for i in D], ms)
        else:
            if ms < len(D) - W:
                nO += 1
                if firstO is None:
                    firstO = (f, [S[i] for i in D], ms)
    verdict = "EXCLUDED" if (nE or nO) else "passes"
    if verbose:
        print(f"{label} m={m} W={W} M={M}: tight {len(tight)}, even-test failures {nE}, odd-test failures {nO} -> {verdict}"
              + (f"; first even failure: support dirs {firstE[1]} (|D|={len(firstE[1])}, matching {firstE[2]})" if firstE else "")
              + (f"; first odd failure: {firstO}" if firstO else ""), flush=True)
    return (verdict, nE, nO)


SURVIVORS = {
    (4, 12): [
        [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15],
    ],
    (4, 13): [
        [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 15, 15, 15],
        [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15],
    ],
}

CALIB = [
    ("plain 1in3 b3", 3, [1, 2, 4], 1),
    ("plain 2in6 b3", 3, [1, 1, 2, 2, 4, 4], 2),
    ("plain 3in9 b3", 3, [1, 2, 4] * 3, 3),
    ("plain 2in8 b4", 4, [1, 2, 4, 8] * 2, 2),
    ("plain 4in16 b4", 4, [1, 2, 4, 8] * 4, 4),
    ("plain 3in15 b5", 5, [1, 2, 4, 8, 16] * 3, 3),
    ("1in4 (7,11,13,14) b4", 4, [7, 11, 13, 14], 1),
    ("a7 2in8 two reps b4", 4, [1, 7, 2, 11, 4, 13, 8, 14], 2),
    ("a8 4in16 2x2 reps b4", 4, [1, 7, 2, 11, 4, 13, 8, 14] * 2, 4),
    ("U3 2in7 b3", 3, list(range(1, 8)), 2),
    ("U3^(2) 4in14 b3", 3, list(range(1, 8)) * 2, 4),
    ("U4 4in15 b4 [excluded by Thm 1]", 4, list(range(1, 16)), 4),
    ("U4 + e1: 4in16 b4 [excluded by Thm 1]", 4, list(range(1, 16)) + [1], 4),
    ("U4 + e1 + e2: 4in17 b4 [excluded by Thm 1]", 4, list(range(1, 16)) + [1, 2], 4),
    ("2in5 b3 [impossible]", 3, [1, 2, 4, 3, 7], 2),
]

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("calib", "all"):
        for lab, m, S, W in CALIB:
            c = clique_test(m, S, W, verbose=False)
            e = epattern_test(m, S, W, verbose=True, label=lab + f" [clique {c}]")
            if e[0] == "deletion":
                print(f"{lab}: deletion test fails ({e[1]} < W)")
    if what in ("surv", "all"):
        for (m, W), L in SURVIVORS.items():
            for S in L:
                c = clique_test(m, S, W, verbose=False)
                extra = [z for z in range(1, 16) if S.count(z) > 3]
                epattern_test(m, S, W, verbose=True, label=f"survivor dirs with >3 copies {extra} [clique {c}]")
