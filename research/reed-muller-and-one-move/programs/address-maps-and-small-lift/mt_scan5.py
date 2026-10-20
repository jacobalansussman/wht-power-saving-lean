"""mt_scan5.py (rm-theory): block 5, the free-target U_5 (8, 9, 10 arrays in 31 moves, one per direction), canonical
Singer targets T_k = sigma^(k) sigma^(k)T (k in Z_31). Exact tests only (deletion, tight forms, tight cliques), on
random target sets with distinct k, plus the all-I slice control (must be excluded). Prints the sets with no tight form.
"""
import sys, random
from mt_lm import Setup
W = int(sys.argv[1]); N = int(sys.argv[2]); seed = int(sys.argv[3])
rng = random.Random(seed)
st = Setup(5, "I:%d" % W, "all", "lex"); print("control all-I:", st.exact_tests(), flush=True)
notight = []; stats = {}
for it in range(N):
    ks = sorted(rng.sample(range(31), W - 1)); ks = [0] + [k for k in ks if k != 0][:W - 1]
    while len(ks) < W:
        k = rng.randrange(1, 31)
        if k not in ks: ks.append(k)
    st = Setup(5, "k:" + ",".join(map(str, sorted(ks))), "all", "lex")
    et = st.exact_tests()
    key = (et["excluded"], et["tight"] > 0)
    stats[key] = stats.get(key, 0) + 1
    if et["tight"] == 0 and len(notight) < 5: notight.append(sorted(ks))
print("W", W, "sets", N, "(excluded, has tight form) -> count:", stats)
print("examples with no tight form:", notight)
