"""vcert.py (v-r3-3in14) [argument local: certificates whose tight pairs lie inside D_F, for the Groebner step]: for each of MY 44 sets, the zero-block certificate with the SMALLEST |D_F| (my search
over all 2^15 forms), re-checked pair by pair; then the list of all forms with support inside U = D_F.
Output: vcerts.json, out_vcert.txt."""
import json
import numpy as np
from vcore import Inst, dir_vec, tgt_vec, EXP, max_clique
import sys
LOCAL = len(sys.argv) > 1 and sys.argv[1] == 'local'   # local: every pair of Y in a tight support INSIDE D_F
sets = json.load(open("vsets.json")); out = []; res = []
for n, s in enumerate(sets):
    I = Inst([tgt_vec(EXP[(2 * k) % 31]) for k in s["ks"]], [dir_vec(x) for x in s["moves"]])
    ts = I.tight_supports(); M = I.M
    J = [[any(t >> i & 1 and t >> j & 1 for t in ts) for j in range(M)] for i in range(M)]
    live = [i for i in range(M) if J[i][i]]
    best = None
    for q in np.nonzero((I.cnt + I.r) < 6)[0]:
        d = int(I.Dmask[q]); Y = max_clique([i for i in live if d >> i & 1], J)
        inside = all(any(t >> i & 1 and t >> j & 1 and t & ~d == 0 for t in ts) for i in Y for j in Y)
        if int(I.cnt[q]) + int(I.r[q]) < 2 * len(Y) and (inside or not LOCAL):
            if best is None or int(I.cnt[q]) < best[1]: best = (int(q), int(I.cnt[q]), int(I.r[q]), Y, d)
    q, dsz, rF, Y, d = best
    # re-check: E vanishes on Y x Y (each pair inside one tight support), F's support contains Y
    assert all(J[i][j] for i in Y for j in Y) and all(d >> i & 1 for i in Y) and dsz + rF < 2 * len(Y)
    U = [i for i in range(M) if d >> i & 1]
    sups = {}
    for qq in range(1, 1 << 15):
        dm = int(I.Dmask[qq])
        if dm and dm & ~d == 0:
            p = tuple(int(I.Pmask[qq]) >> l & 1 for l in range(3))
            assert sups.get(dm, p) == p
            sups[dm] = p
    supl = [([i for i in U if dm >> i & 1], p) for dm, p in sorted(sups.items())]
    tight_in = [t for t in ts if t & ~d == 0]
    res.append({"n": n, "ks": s["ks"], "moves": s["moves"], "U": U, "U_dirs": [s["moves"][i] for i in U],
                "Y": Y, "rF": rF, "sups": supl})
    s["U_dirs"] = [s["moves"][i] for i in U]
    out.append(f"set {n} ks={tuple(s['ks'])}: |D_F|={dsz} r_F={rF} Y={Y} U={U} dirs={[s['moves'][i] for i in U]} "
               f"forms-with-support-in-U: {len(supl)} distinct supports; tight supports inside U: {len(tight_in)}")
from collections import Counter
out.append("sizes |D_F|: " + str(Counter(len(r["U"]) for r in res)) + "; number of supports in U: " +
           str(Counter(len(r["sups"]) for r in res)) + "; repeated direction inside U: " +
           str(Counter(len(set(r["U_dirs"])) < len(r["U_dirs"]) for r in res)))
json.dump(res, open("vcerts_local.json" if LOCAL else "vcerts.json", "w")); json.dump(sets, open("vsets.json", "w"))
open("out_vcert_local.txt" if LOCAL else "out_vcert.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out[-3:]))
