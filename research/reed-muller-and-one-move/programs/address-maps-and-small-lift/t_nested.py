"""t_nested.py (rm-theory): the NESTED-SUPPORT test for programs given by a direction list and per-array targets:
for forms F, G with D_F a subset of D_G (as multisets of moves),  |D_F| + r_F <= |D_G| + r_G.
Proof: rank E[D,D] = |D| - r (NOTE 8.5) and rank of a matrix <= rank of a principal block + 2 x (rows outside it).
Usage: python3 t_nested.py m Wspec dirs   (mt_lm.Setup syntax: k:..., gl:..., I:n ; all | frames | list:...)
       python3 t_nested.py scan m W N seed   (random distinct canonical targets, one move per direction)"""
import sys, random
import numpy as np
from mt_lm import Setup


def nested(st, W):
    vals = st.vals.astype(np.int64); M = vals.shape[1]
    cnt = vals.sum(axis=1); r = st.phi.sum(axis=1)
    # multiset supports: move j in D_F; as bitmask over the M moves (M <= 62)
    masks = (vals * (1 << np.arange(M, dtype=np.int64))[None, :]).sum(axis=1)
    key = cnt + r
    order = np.argsort(cnt)
    # only pairs with |D_G| - |D_F| < r_F - r_G <= W can fail
    best = None
    um, idx = np.unique(masks, return_index=True)        # distinct supports (r is a function of the support)
    cntu = cnt[idx]; ru = r[idx]; keyu = cntu + ru
    nz = cntu > 0
    um, cntu, ru, keyu = um[nz], cntu[nz], ru[nz], keyu[nz]
    for i in np.argsort(cntu):
        # candidates G: superset of F, with key_G < key_F, and |D_G| < |D_F| + W
        sel = (cntu > cntu[i]) & (cntu < cntu[i] + W + 1) & (keyu < keyu[i])
        if not sel.any(): continue
        sup = (um[sel] & um[i]) == um[i]
        if sup.any():
            j = np.where(sel)[0][np.where(sup)[0][0]]
            return dict(excluded=True, F=(int(cntu[i]), int(ru[i])), G=(int(cntu[j]), int(ru[j])))
    return dict(excluded=False)


if __name__ == "__main__":
    if sys.argv[1] == "scan":
        m, W, N, seed = map(int, sys.argv[2:6]); n = (1 << m) - 1
        rng = random.Random(seed); stats = {}
        for it in range(N):
            ks = [0] + sorted(rng.sample(range(1, n), W - 1))
            st = Setup(m, "k:" + ",".join(map(str, ks)), "all", "lex")
            et = st.exact_tests(); nt = nested(st, W)
            key = ("old-excluded" if et["excluded"] else ("tight" if et["tight"] else "no-tight"),
                   "NESTED-excl" if nt["excluded"] else "nested-pass")
            stats[key] = stats.get(key, 0) + 1
            if nt["excluded"] and not et["excluded"] and stats[key] <= 3:
                print("   excluded by nested:", ks, nt, flush=True)
        print(f"m={m} W={W} in {n}: {N} random canonical target sets (seed {seed}):", stats, flush=True)
    else:
        m = int(sys.argv[1]); st = Setup(m, sys.argv[2], sys.argv[3], "lex")
        print(sys.argv[1:], st.exact_tests(), nested(st, st.W), flush=True)
