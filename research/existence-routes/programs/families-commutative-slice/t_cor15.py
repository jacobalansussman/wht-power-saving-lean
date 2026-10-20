# check of Corollary 15's combinatorics with the quadratic-form tables (independent of Lemma 12's hand proof):
# F = {z1 = z2 = 1}; h = z3, h2 = z4 on F; each of the six sets H, H', K u K', K u L', L u K', L u L' must lie inside
# the support of an ODD binary quadratic form of minimum weight 2^(m-2) (an odd (m-2)-flat avoiding 0); F itself must
# be the support of an EVEN form.
import sys, numpy as np
from refl import quad_forms
for m in (4, 5, 6):
    vals, tr = quad_forms(m)
    bit = lambda z, j: (z >> j) & 1
    F = [z for z in range(1, 1 << m) if bit(z, 0) and bit(z, 1)]
    h = lambda z: bit(z, 2); h2 = (lambda z: bit(z, 3)) if m >= 5 else (lambda z: bit(z, 2) ^ 1 if False else bit(z, 3) if m > 3 else 0)
    if m == 4:  # F is a plane {z1=z2=1}: coordinates z3, z4 free; quarters are single points
        h = lambda z: bit(z, 2); h2 = lambda z: bit(z, 3)
    sets = {"H": [z for z in F if h(z) == 0], "H'": [z for z in F if h(z) == 1],
            "KuK'": [z for z in F if h2(z) == 0], "KuL'": [z for z in F if h(z) ^ h2(z) == 0],
            "LuK'": [z for z in F if h(z) ^ h2(z) == 1], "LuL'": [z for z in F if h2(z) == 1]}
    wt = vals[:, 1:].sum(1)
    Ff = [f for f in range(len(vals)) if set(np.nonzero(vals[f])[0].tolist()) == set(F)]
    out = [f"m={m}: F = {F} is the support of form(s) with trace {[int(tr[f]) for f in Ff]} (0 = even)"]
    for name, T in sets.items():
        ok = [f for f in range(len(vals)) if tr[f] == 1 and wt[f] == 2 ** (m - 2) and all(vals[f, z] for z in T)]
        out.append(f"  {name} = {T}: odd minimum-weight forms containing it: {len(ok)}")
    print("\n".join(out), flush=True)
