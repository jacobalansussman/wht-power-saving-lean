"""tr_check.py: re-test a survivors file 'dirs|targets' with the leaf tests of t_small.py (TERMRANK=1 adds the term-rank
test). Prints how many survive and why the others fail."""
import sys, os
import numpy as np
from t_free import forms
from t_small import targets_of, leaf_tests
m = int(sys.argv[1]); fn = sys.argv[2]
pairs, mon, coeffs, vals = forms(m)
Tm = targets_of(m, pairs)
par = (coeffs @ Tm.T % 2).T.astype(np.int16)
ok = 0; bad = {}
for line in open(fn):
    d, t = line.strip().split("|")
    ms = [int(x) for x in d.split(",")]; ts = [int(x) for x in t.split(",")]
    cnt = vals[:, ms].sum(axis=1).astype(np.int16)
    r = par[ts].sum(axis=0)
    assert (cnt >= r).all()
    res = leaf_tests(vals[:, ms].astype(bool), cnt, r, len(ts))
    if res[0]: ok += 1
    else:
        key = "clique" if res[3] is None else ("zero-block" if res[3][3] >= 0 else "term-rank")
        bad[key] = bad.get(key, 0) + 1
        if bad[key] <= 2: print("  excluded", key, line.strip(), res[3])
print(f"{fn}: {ok} survive, excluded {bad}  (TERMRANK={os.environ.get('TERMRANK','0')})")
