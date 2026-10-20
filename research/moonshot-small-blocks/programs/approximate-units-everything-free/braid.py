"""braid.py -- xw-border.  The braid relation D_z S_z D_z = S_z D_z S_z at exact Q (S_z = i^(z.x) per column, free in
XW), and whether a fibrewise sub-call keeps it (it should not)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from keptlib import Fib
rng = np.random.default_rng(11)
for f in range(1, 5):
    F = Fib(f)
    v = rng.standard_normal((F.n, 2)) + 1j * rng.standard_normal((F.n, 2))
    Y = F.rand_Y(rng)
    out = []
    for z, S in [(1, "S1"), (2, "S2")]:
        lhs = F.D(z, F.free(S, F.D(z, v))); rhs = F.free(S, F.D(z, F.free(S, v)))
        lg = F.D(z, F.free(S, F.D(z, v, Y)), Y); rg = F.free(S, F.D(z, F.free(S, v), Y))
        out.append("z=%d exact %.1e generic %.2f" % (z, np.linalg.norm(lhs - rhs) / np.linalg.norm(v), np.linalg.norm(lg - rg) / np.linalg.norm(v)))
    print("f=%d: " % f + " | ".join(out))
