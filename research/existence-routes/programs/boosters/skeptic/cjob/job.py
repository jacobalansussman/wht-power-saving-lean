# v-e3-boost Colab job: (1) booster LP at registers 5 and 6 with my certificates on M = m .. m+6 (spare rows), with and
# without the quoted hooks; (2) S6, S7, S8 against my families on M = m .. 14.
import sys, time
import numpy as np
import vsurv as S, vlp
t0 = time.time()
named = {"S6": {2: -1, 5: -4, 3: 2, 4: 1, 6: 2}, "S7": {2: -1, 5: -1, 6: -1, 3: 2, 7: 1},
         "S8": {2: -1, 5: -1, 6: -2, 3: 1, 4: 2, 8: 1}}
for name, v in named.items():
    m = max(v); res = S.tests(v, m, Mextra=14 - m); bad = [r for r in res if r[1] < 0]
    print(f"{name}: {len(res)} tests on M = {m}..14 (+ b=2), failed {len(bad)} {bad[:4]}  [{time.time()-t0:.0f}s]", flush=True)
grid = [i / 50 for i in range(1, 50)] + [0.9983]
for m, Mx in ((5, 6), (6, 6)):
    for hooks in (False, True):
        R = vlp.rows(m, hooks, Mextra=Mx)
        neg = []; worst = (0, None, None)
        for t in grid:
            f, vv = vlp.boost_min(m, t, R)
            if f < -1e-9:
                neg.append(t)
                if f < worst[0]: worst = (f, t, vv)
        print(f"m={m} M<= {m+Mx} hooks={hooks}: {len(R)} rows; booster-negative at {len(neg)}/{len(grid)}"
              + (f" ({min(neg)}..{max(neg)}); worst {worst[0]:.3g} at {worst[1]}, v = {np.round(worst[2] / np.abs(worst[2]).max(), 4).tolist()}" if neg else "")
              + f"  [{time.time()-t0:.0f}s]", flush=True)
