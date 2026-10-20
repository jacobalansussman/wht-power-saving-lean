"""orbit_count.py (rm2-phase): sum of orbit sizes (under the direction multiset's stabiliser) of the survivor classes
in a fam_enum dump, split by the number of alternating targets; cross-check against t_small_alt.py tuple counts."""
import sys
import numpy as np
from t_free import gl_elements, apply as gl_apply
from tile_alt import target_matrices, mat_key, perm_table
m = int(sys.argv[1]); path = sys.argv[2]
mats, nonalt = target_matrices(m); keyidx = {mat_key(S): i for i, S in enumerate(mats)}
rows = [L.strip().split("|") for L in open(path) if "|" in L]
ms = tuple(int(z) for z in rows[0][0].split(","))
Gs = [h for h in gl_elements(m) if tuple(sorted(gl_apply(h, z) for z in ms)) == ms]
perm = perm_table(Gs, mats, keyidx, m)
tot = {}
for _, t in rows:
    tt = np.array([int(x) for x in t.split(",")])
    imgs = set(tuple(sorted(perm[g, tt])) for g in range(len(Gs)))
    na = int((tt >= nonalt).sum())
    tot[na] = tot.get(na, 0) + len(imgs)
print(f"dirs {ms}: stabiliser {len(Gs)}; {len(rows)} classes; tuples by number of alternating targets {dict(sorted(tot.items()))}; "
      f"with at least one alternating: {sum(v for k, v in tot.items() if k > 0)}; total {sum(tot.values())}")
