"""mt_scan5z.py (rm-theory, after the restart): mt_scan5.py's random canonical-Singer target sets for the uniform family
at block 5 (W in 31, one move per direction), now with the zero-block test of NOTE 8.7 added (t_small.leaf_tests)."""
import sys, random
import numpy as np
from mt_lm import Setup
from t_small import leaf_tests
W = int(sys.argv[1]); N = int(sys.argv[2]); seed = int(sys.argv[3])
rng = random.Random(seed); stats = {}
for it in range(N):
    ks = [0] + sorted(rng.sample(range(1, 31), W - 1))
    st = Setup(5, "k:" + ",".join(map(str, ks)), "all", "lex")
    et = st.exact_tests()
    cnt = st.vals.sum(axis=1).astype(np.int16); r = st.phi.sum(axis=1).astype(np.int16)
    if et["excluded"]:
        key = "excluded by deletion/cliques"
    elif et["tight"] == 0:
        key = "no tight form"
    else:
        res = leaf_tests(st.vals.astype(bool), cnt, r, W)
        key = "tight forms, pass Z" if res[0] else "tight forms, EXCLUDED by Z"
        if not res[0]: print("  Z excludes k =", ks, res[3], flush=True)
    stats[key] = stats.get(key, 0) + 1
print(f"W={W} in 31, {N} random canonical target sets (seed {seed}):", stats)
