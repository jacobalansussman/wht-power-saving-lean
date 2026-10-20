"""k2_sumzero_nonalt.py (skeptic-2): the only 4-in-15 "cores" (each of the 15 directions of F_2^4 once, 4 arrays)
that L-det allows WITHOUT scratch are those whose four targets sum to 0. k1 showed that none of rm2-phase A3(d)'s
1,612 survivors with an alternating target sums to 0. Here: every 4-multiset of NON-alternating nonsingular targets
with sum 0, through rm2-phase's own tests S+D+C+Z (t_small_alt.leaf_tests, imported unchanged; TERMRANK as there).

GL(4,2) fixes the direction set and is transitive on the 420 non-alternating targets, so T1 = I loses nothing.
"""
import sys, os, time
import numpy as np
RM = os.environ.get('RM_PHASE_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                   'reed-muller-and-one-move', 'programs', 'phase-letters', 'ends-only-block4'))
RMD = os.environ.get('RM_PHASE_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                   'reed-muller-and-one-move', 'data', 'phase-letters', 'ends-only-block4'))
sys.path.insert(0, RM)
import t_small_alt
from t_small_alt import leaf_tests
from t_free import forms as tf_forms
from tile_alt import target_matrices, mat_key

t0 = time.time()
m, W = 4, 4
t_small_alt.ALT = True
pairs, mon, coeffs, vals = tf_forms(m)
Tm = t_small_alt.targets_of(m, pairs); nonalt = t_small_alt.NONALT
mats, _ = target_matrices(m)
key = {mat_key(S): i for i, S in enumerate(mats)}
par = (coeffs @ Tm.T % 2).T.astype(np.int16)          # targets x forms
ms = list(range(1, 16))
cnt = np.zeros(vals.shape[0], dtype=np.int16)
for z in ms: cnt += vals[:, z].astype(np.int16)
Zv = vals[:, ms].astype(bool)
print("forms:", vals.shape[0], " support sizes:", sorted(set(int(c) for c in cnt)), " TERMRANK =", t_small_alt.TERMRANK)
I = [[int(i == j) for j in range(m)] for i in range(m)]
i0 = key[mat_key(I)]
def add(A, B): return [[a ^ b for a, b in zip(r, s)] for r, s in zip(A, B)]

# calibration: a plain 4 in 16 (every target I) must pass D on 16 moves; 5 arrays on the 15 directions must fail.
seen = set(); stats = {"tuples": 0, "fail_D": 0, "fail_CZ": 0, "pass": 0}
survivors = []
for a in range(nonalt):
    for b in range(a, nonalt):
        S4 = add(add(I, mats[a]), mats[b])
        k4 = key.get(mat_key(S4))
        if k4 is None or k4 >= nonalt: continue        # T4 must be nonsingular and non-alternating
        tup = tuple(sorted((i0, a, b, k4)))
        if tup in seen: continue
        seen.add(tup); stats["tuples"] += 1
        r = par[list(tup)].sum(axis=0)
        assert (r % 2 == 0).all()                       # sum zero -> every form is even
        if not (cnt >= r).all(): stats["fail_D"] += 1; continue
        res = leaf_tests(Zv, cnt, r, W)
        if res[0]:
            stats["pass"] += 1; survivors.append((tup, res[1], res[2]))
        else:
            stats["fail_CZ"] += 1
print("sum-zero non-alternating 4-multisets containing I:", stats, f"({time.time()-t0:.1f} s)")
nt0 = sum(1 for s in survivors if s[1] == 0)
print("survivors with no tight form:", nt0, " with tight forms:", len(survivors) - nt0)
for s in survivors[:12]: print("  ", s)
