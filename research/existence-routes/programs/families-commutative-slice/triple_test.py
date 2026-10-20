"""triple_test.py (ex2-families, session 2, 14 Oct 2026 23:10 EDT): the slack-one triple test (PROOF.md Theorem 12).

E = I + Gamma - K.  Facts: even Q -> E[D_Q, D_Q] nonsingular; odd Q -> rank E[D_Q, D_Q] = |D_Q| - W (its slack).
Lemma 11: an n x n matrix (n >= 2) whose principal 2x2 and 3x3 submatrices all have rank <= 1 is singular.
TEST: if some even Q with |D_Q| >= 2 has every set of <= 3 of its moves inside the support of an odd form of slack <= 1,
then no commutative-slice program exists.  Copies of a direction lie in exactly the same supports, so the test works on
direction sets (with multiplicity: a set of 3 moves of D_Q uses at most as many copies of a direction as D_Q holds).
usage: triple_test.py m W dirs...   or   triple_test.py calib | family
"""
import sys, itertools
import numpy as np
from refl import quad_forms

_cache = {}


def forms(m):
    if m not in _cache:
        _cache[m] = quad_forms(m)
    return _cache[m]


def triple_test(m, S, W, verbose=False, label="", max_even=None):
    vals, tr = forms(m)
    S = list(S)
    dirs = sorted(set(S))
    mult = {z: S.count(z) for z in dirs}
    cnt = np.array([mult.get(z, 0) for z in range(1 << m)])
    wt = (vals * cnt[None, :]).sum(axis=1)         # |D_Q| with multiplicity
    odd = tr == 1
    if wt[odd].min() < W:
        return ("deletion", int(wt[odd].min()))
    slack1 = np.nonzero(odd & (wt <= W + 1))[0]
    covered = set()
    for f in slack1:
        D = [z for z in dirs if vals[f, z]]
        for k in (1, 2, 3):
            for T in itertools.combinations(D, k):
                covered.add(T)
    even = np.nonzero((~odd) & (wt >= 2))[0]
    even = even[np.argsort(wt[even], kind="stable")]
    if max_even:
        even = even[:max_even]
    hit = None
    for f in even:
        D = [z for z in dirs if vals[f, z]]
        ok = True
        # sets of <= 3 moves of D: direction sets T of size 1..3 (size-1 / size-2 sets need enough copies only when the
        # move set has repeated directions; covering T as a direction set covers every move set with those directions)
        for k in (1, 2, 3):
            for T in itertools.combinations(D, k):
                if T not in covered:
                    ok = False; break
            if not ok:
                break
        if ok:
            hit = (int(f), D, int(wt[f])); break
    verdict = "EXCLUDED" if hit else "passes"
    if verbose:
        print(f"{label} m={m} W={W} M={len(S)}: odd forms of slack <= 1: {len(slack1)}; even forms checked {len(even)}"
              f" -> {verdict}" + (f"; even support fully triple-covered: dirs {hit[1]} (|D| = {hit[2]})" if hit else ""),
              flush=True)
    return (verdict, hit)


CALIB = [
    ("plain 1in3 b3", 3, [1, 2, 4], 1),
    ("plain 2in6 b3", 3, [1, 1, 2, 2, 4, 4], 2),
    ("plain 3in9 b3", 3, [1, 2, 4] * 3, 3),
    ("plain 2in8 b4", 4, [1, 2, 4, 8] * 2, 2),
    ("plain 4in16 b4", 4, [1, 2, 4, 8] * 4, 4),
    ("plain 3in15 b5", 5, [1, 2, 4, 8, 16] * 3, 3),
    ("plain 1in5 b5", 5, [1, 2, 4, 8, 16], 1),
    ("1in4 (7,11,13,14) b4", 4, [7, 11, 13, 14], 1),
    ("a7 2in8 two reps b4", 4, [1, 7, 2, 11, 4, 13, 8, 14], 2),
    ("a8 4in16 2x2 reps b4", 4, [1, 7, 2, 11, 4, 13, 8, 14] * 2, 4),
    ("U3 2in7 b3", 3, list(range(1, 8)), 2),
    ("U3 1in7 b3 (W below)", 3, list(range(1, 8)), 1),
    ("U3^(2) 4in14 b3", 3, list(range(1, 8)) * 2, 4),
    ("U3^(2) 3in14 b3", 3, list(range(1, 8)) * 2, 3),
    ("U4 3in15 b4 (exists? 3 plain arrays fit in 12 of its moves)", 4, list(range(1, 16)), 3),
    ("U4 4in15 b4 [excluded by Thm 1]", 4, list(range(1, 16)), 4),
]

if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "calib":
        for lab, m, S, W in CALIB:
            triple_test(m, S, W, verbose=True, label=lab)
    elif a and a[0] == "family":
        # U_m^(r) with W = r 2^(m-2) - s
        for (m, r, s) in [(4, 1, 1), (4, 2, 1), (4, 3, 1), (4, 5, 1), (4, 8, 1), (5, 1, 1), (5, 1, 2), (5, 2, 1), (6, 1, 1),
                          (6, 1, 2), (6, 1, 3)]:
            W = r * 2 ** (m - 2) - s
            S = list(range(1, 1 << m)) * r
            M = len(S)
            sv = np.log(m * W / M) / np.log(m)
            triple_test(m, S, W, verbose=True, label=f"U_{m}^({r}) W={W} (slack {s}; saving if it existed {sv:.3f})")
    else:
        m, W = int(a[0]), int(a[1]); S = list(map(int, a[2:]))
        triple_test(m, S, W, verbose=True)
