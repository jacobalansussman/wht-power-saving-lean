"""cover_test.py (ex2-families, session 2): the slack-s covering test (PROOF.md Theorem 12, general s).

Lemma 11 (general): if every principal submatrix of order <= s+2 of an n x n matrix (n >= s+1) has rank <= s, the
matrix is singular.  With E = I + Gamma - K (even Q: E[D_Q, D_Q] nonsingular; odd Q: rank E[D_Q, D_Q] = slack(Q)):
TEST(s): if some even Q with |D_Q| >= s+1 has every set of <= s+2 of its moves inside the support of an odd form of
slack <= s, then no commutative-slice program exists.  (s = 1 is triple_test.py.)
usage: cover_test.py s m W dirs...    or    cover_test.py family [maxeven]
"""
import sys, itertools
import numpy as np
from refl import quad_forms

_cache = {}


def forms(m):
    if m not in _cache:
        _cache[m] = quad_forms(m)
    return _cache[m]


def cover_test(m, S, W, s, verbose=False, label="", max_even=None, max_even_weight=None):
    vals, tr = forms(m)
    S = list(S)
    dirs = sorted(set(S))
    cnt = np.zeros(1 << m, dtype=int)
    for z in S:
        cnt[z] += 1
    wt = (vals * cnt[None, :]).sum(axis=1)
    odd = tr == 1
    if wt[odd].min() < W:
        return ("deletion", int(wt[odd].min()))
    good = np.nonzero(odd & (wt <= W + s))[0]
    k = s + 2
    covered = set()
    for f in good:
        D = [z for z in dirs if vals[f, z]]
        for j in range(1, min(k, len(D)) + 1):
            covered.update(itertools.combinations(D, j))
    even = np.nonzero((~odd) & (wt >= s + 1))[0]
    if max_even_weight:
        even = even[wt[even] <= max_even_weight]
    even = even[np.argsort(wt[even], kind="stable")]
    if max_even:
        even = even[:max_even]
    hit = None
    for f in even:
        D = [z for z in dirs if vals[f, z]]
        ok = True
        for j in range(1, min(k, len(D)) + 1):
            for T in itertools.combinations(D, j):
                if T not in covered:
                    ok = False; break
            if not ok:
                break
        if ok:
            hit = (int(f), D, int(wt[f])); break
    verdict = "EXCLUDED" if hit else "passes"
    if verbose:
        print(f"{label} m={m} W={W} M={len(S)} s={s}: odd forms of slack <= {s}: {len(good)}; even forms checked "
              f"{len(even)} -> {verdict}" + (f"; covered even support: dirs {hit[1]} (|D| = {hit[2]})" if hit else ""),
              flush=True)
    return (verdict, hit)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "family":
        mx = int(a[1]) if len(a) > 1 else None
        for (m, sl) in [(4, [1, 2]), (5, [1, 2, 3]), (6, [1, 2, 3, 4]), (7, [1, 2, 3])]:
            for s in sl:
                W = 2 ** (m - 2) - s
                S = list(range(1, 1 << m))
                sv = np.log(m * W / len(S)) / np.log(m)
                cover_test(m, S, W, s, verbose=True, max_even=mx, max_even_weight=2 ** (m - 2),
                           label=f"U_{m} W={W} (slack {s}; saving if it existed {sv:.3f})")
    else:
        s, m, W = int(a[0]), int(a[1]), int(a[2]); S = list(map(int, a[3:]))
        cover_test(m, S, W, s, verbose=True)
