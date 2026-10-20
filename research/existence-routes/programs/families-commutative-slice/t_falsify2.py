# falsification test for Theorems 10 and 12 (even test, odd term-rank test, slack-one triple test): multisets that
# certainly carry a unit must never be EXCLUDED.  Units: disjoint unions of W representations of I_m (one register per
# array; plain type), optionally with extra unused moves added (adding an idle move keeps a unit a unit).
# usage: t_falsify2.py m ntests seed extra_max
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
from even_test import epattern_test
from triple_test import triple_test
from clique_test import clique_test
m, ntests, seed, extra_max = map(int, sys.argv[1:5])
rng = np.random.default_rng(seed)
def cls(z):
    v = [(z >> j) & 1 for j in range(m)]
    return tuple(v[j] * v[k] for j in range(m) for k in range(j, m))
I = tuple(1 if j == k else 0 for j in range(m) for k in range(j, m))
pool = []
while len(pool) < 300:
    k = int(rng.integers(m, m + 5))
    zs = [int(z) for z in rng.integers(1, 1 << m, size=k)]
    s = np.zeros(len(I), dtype=int)
    for z in zs: s = (s + np.array(cls(z))) % 2
    if tuple(s.tolist()) == I: pool.append(zs)
bad = 0; nt = ne = 0
for t in range(ntests):
    W = int(rng.integers(1, 5))
    S = sum((pool[i] for i in rng.integers(0, len(pool), size=W)), [])
    S += [int(z) for z in rng.integers(1, 1 << m, size=int(rng.integers(0, extra_max + 1)))]
    rng.shuffle(S)
    c = clique_test(m, S, W, verbose=False)
    e = epattern_test(m, S, W)
    tt = triple_test(m, S, W)
    if c[0] == "EXCLUDED" or e[0] == "EXCLUDED" or tt[0] == "EXCLUDED":
        bad += 1; print("FALSIFIED?", m, W, S, c, e, tt, flush=True)
    # count how often the tests had something to work with
    vals_has_tight = c[1] >= 1 if c[0] == "passes" else False
    nt += vals_has_tight
print(f"m={m} extra<= {extra_max}: {ntests} existing units tested, EXCLUDED by clique/even/odd/triple: {bad}; "
      f"cases with at least one tight form: {nt}", flush=True)
