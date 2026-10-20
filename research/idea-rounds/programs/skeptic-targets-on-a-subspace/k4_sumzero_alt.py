"""k4_sumzero_alt.py (skeptic-2): independent re-check of k1's consequence, without A3(d)'s survivor file.
Every 4-multiset of nonsingular symmetric 4x4 targets with AT LEAST ONE alternating target and sum 0, on the 15
directions of F_2^4 each used once, through D + rm2-phase's leaf_tests (C, Z).
GL(4,2) fixes the direction set and is transitive on the 28 alternating targets (20160/720), so T1 = one fixed
alternating target loses nothing.
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
t0 = time.time(); m, W = 4, 4
t_small_alt.ALT = True
pairs, mon, coeffs, vals = tf_forms(m)
Tm = t_small_alt.targets_of(m, pairs); nonalt = t_small_alt.NONALT
mats, _ = target_matrices(m); key = {mat_key(S): i for i, S in enumerate(mats)}
par = (coeffs @ Tm.T % 2).T.astype(np.int16)
ms = list(range(1, 16))
cnt = np.zeros(vals.shape[0], dtype=np.int16)
for z in ms: cnt += vals[:, z].astype(np.int16)
Zv = vals[:, ms].astype(bool)
def add(A, B): return [[a ^ b for a, b in zip(r, s)] for r, s in zip(A, B)]
a0 = nonalt                                            # one alternating target
seen = set(); stats = {"tuples": 0, "by_nalt": {}, "fail_D": 0, "fail_CZ": 0, "pass": 0}
for a in range(len(mats)):
    for b in range(a, len(mats)):
        S4 = add(add(mats[a0], mats[a]), mats[b])
        k4 = key.get(mat_key(S4))
        if k4 is None: continue                         # T4 must be nonsingular
        tup = tuple(sorted((a0, a, b, k4)))
        if tup in seen: continue
        seen.add(tup); stats["tuples"] += 1
        na = sum(1 for x in tup if x >= nonalt); stats["by_nalt"][na] = stats["by_nalt"].get(na, 0) + 1
        r = par[list(tup)].sum(axis=0)
        if not (cnt >= r).all(): stats["fail_D"] += 1; continue
        res = leaf_tests(Zv, cnt, r, W)
        if res[0]: stats["pass"] += 1; print("  PASS", tup, res[1], res[2])
        else: stats["fail_CZ"] += 1
print("sum-zero 4-multisets with >= 1 alternating target (one alternating fixed):", stats, f"({time.time()-t0:.1f} s)")
