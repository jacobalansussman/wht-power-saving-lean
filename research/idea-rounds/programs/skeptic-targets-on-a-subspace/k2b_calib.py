"""k2b_calib.py (skeptic-2): calibrate k2's harness.
(1) All 1,612 rm2-phase A3(d) survivors must PASS the same D + leaf_tests call (they passed fam_enum).
(2) Plain controls must pass: 4 in 16 = four arrays each owed I on the standard basis (16 moves, e1..e4 four times);
    and 3 in 15 with sum-zero targets I, S, I+S plus ... (skipped).  Negative: 5 in 15 with five copies of I must fail D/CZ.
"""
import sys, os
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
m = 4
t_small_alt.ALT = True
pairs, mon, coeffs, vals = tf_forms(m)
Tm = t_small_alt.targets_of(m, pairs); nonalt = t_small_alt.NONALT
mats, _ = target_matrices(m); key = {mat_key(S): i for i, S in enumerate(mats)}
par = (coeffs @ Tm.T % 2).T.astype(np.int16)
def test(ms, tup, W):
    cnt = np.zeros(vals.shape[0], dtype=np.int16)
    for z in ms: cnt += vals[:, z].astype(np.int16)
    Zv = vals[:, list(ms)].astype(bool)
    r = par[list(tup)].sum(axis=0)
    if not (cnt >= r).all(): return "fail D"
    res = leaf_tests(Zv, cnt, r, W)
    return "pass" if res[0] else "fail CZ"
rows = open(os.path.join(RMD, "surv_fe_alt_415_all.txt")).read().split()
out = {}
for line in rows:
    d, t = line.split("|")
    v = test([int(x) for x in d.split(",")], [int(x) for x in t.split(",")], 4)
    out[v] = out.get(v, 0) + 1
print("A3(d) survivors through this harness:", out)
I = [[int(i == j) for j in range(m)] for i in range(m)]; iI = key[mat_key(I)]
std = [1, 2, 4, 8]                       # direction z encoded as an integer (bit i = coordinate i)
print("plain 4 in 16 (e1..e4 four times, targets I x4):", test(std * 4, [iI] * 4, 4))
print("5 in 15, each direction once, targets I x5:", test(list(range(1, 16)), [iI] * 5, 5))
print("4 in 15, each direction once, targets I x4 (sum zero):", test(list(range(1, 16)), [iI] * 4, 4))
