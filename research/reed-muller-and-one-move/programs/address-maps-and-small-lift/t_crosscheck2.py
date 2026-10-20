"""t_crosscheck2.py: for one direction multiset, every target W-multiset that passes deletion (t_small's code) is
re-tested with mt_lm.Setup.exact_tests (independent deletion/tight/clique code); verdicts compared."""
import sys, itertools
import numpy as np
from t_free import forms
from t_small import targets_of, leaf_tests
from mt_lm import Setup
m, W = int(sys.argv[1]), int(sys.argv[2]); ms = [int(x) for x in sys.argv[3].split(",")]
pairs, mon, coeffs, vals = forms(m)
Tm = targets_of(m, pairs); par = (coeffs @ Tm.T % 2).T.astype(np.int16)
cnt = vals[:, ms].sum(axis=1).astype(np.int16)
inV = [int(i) for i in np.where(~(par[:, cnt == 0].astype(bool)).any(axis=1))[0]]
n = agreeD = agreeC = 0; passC = 0
for ts in itertools.combinations_with_replacement(inV, W):
    r = par[list(ts)].sum(axis=0)
    if not (cnt >= r).all(): continue
    n += 1
    et = Setup(m, "gl:" + ",".join(map(str, ts)), "list:" + ",".join(map(str, ms)), "lex").exact_tests()
    res = leaf_tests(vals[:, ms].astype(bool), cnt, r, W)
    mineC = res[0] or res[3] is not None
    agreeD += int(et["violations"] == 0); agreeC += int(mineC == (not et["excluded"])); passC += int(mineC)
print(f"m={m} W={W} dirs {ms}: {n} target tuples pass my deletion; mt_lm finds no violation in {agreeD}; clique verdicts "
      f"agree in {agreeC}; pass cliques: {passC}")
