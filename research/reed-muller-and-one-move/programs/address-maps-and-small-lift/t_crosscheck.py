"""t_crosscheck.py: compare t_small.py's tests (span, deletion, cliques) with mt_lm.py's independent exact_tests on random
small-lift configurations (block m, W arrays, M moves; targets drawn from those in the span so the comparison is not
trivial)."""
import sys, random
import numpy as np
from t_free import forms
from t_small import targets_of, leaf_tests
from mt_lm import Setup
m, W, M, N, seed = map(int, sys.argv[1:6])
rng = random.Random(seed)
pairs, mon, coeffs, vals = forms(m)
Tm = targets_of(m, pairs)
par = (coeffs @ Tm.T % 2).T.astype(np.int16)
agree = 0; tested = 0; kinds = {}
while tested < N:
    ms = sorted(rng.randrange(1, 1 << m) for _ in range(M))
    cnt = vals[:, ms].sum(axis=1).astype(np.int16)
    Z = cnt == 0
    inV = np.where(~(par[:, Z].astype(bool)).any(axis=1))[0]
    if len(inV) == 0: continue
    ts = [int(rng.choice(inV)) for _ in range(W)]
    r = par[ts].sum(axis=0)
    mine_D = bool((cnt >= r).all())
    mine = mine_D and leaf_tests(vals[:, ms].astype(bool), cnt, r, W)[0] if mine_D else False
    # mine without Z: clique only
    st = Setup(m, "gl:" + ",".join(map(str, ts)), "list:" + ",".join(map(str, ms)), "lex")
    et = st.exact_tests()
    theirs = not et["excluded"]
    res = leaf_tests(vals[:, ms].astype(bool), cnt, r, W) if mine_D else (False, 0, 0, None)
    mine_C = mine_D and (res[0] or res[3] is not None)      # passes deletion and cliques (Z may still fail)
    tested += 1; agree += int(mine_C == theirs)
    k = (theirs, mine)
    kinds[k] = kinds.get(k, 0) + 1
print(f"m={m} W={W} M={M}: {N} random configurations; span+deletion+clique verdicts agree with mt_lm in {agree}; "
      f"(mt_lm pass, mine incl. Z pass) counts: {kinds}")
