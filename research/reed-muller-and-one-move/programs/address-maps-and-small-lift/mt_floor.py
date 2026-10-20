"""mt_floor.py: where does the LM residual sit at its floor? Per form: residual norm^2, the form's weight on the
moves, n_odd, slack, and whether the form is linear / quadratic (b and c parts in the canonical placement)."""
import sys
import numpy as np
from mt_lm import Setup, MTSystem
m = int(sys.argv[1]); Wspec = sys.argv[2]; dspec = sys.argv[3]; order = sys.argv[4]; seed = int(sys.argv[5])
st = Setup(m, Wspec, dspec, order); sy = MTSystem(st)
rng = np.random.default_rng(seed)
x0 = (rng.standard_normal(sy.nvar) + 1j * rng.standard_normal(sy.nvar)) / np.sqrt(2); x0[2 * sy.nA:] *= 0.3
x, f, it = sy.solve(x0, iters=int(sys.argv[6]) if len(sys.argv) > 6 else 400)
r = sy.residual(x).reshape(len(sy.Ds), sy.W, sy.W)
per = np.sum(np.abs(r) ** 2, axis=(1, 2))
print(f"{Wspec} {dspec} {order}: residual {f:.4f} after {it} its")
idx = np.argsort(-per)
from collections import Counter
cnt = Counter()
for i in idx[:40]:
    D = sy.Ds[i]; nodd = int(sy.phis[i].sum())
    print(f"  form#{i}: res {per[i]:.4f}  |D|={len(D)} n_odd={nodd} slack={len(D)-nodd} diag {np.round(np.diag(r[i]).real,3)} offmax {np.max(np.abs(r[i]-np.diag(np.diag(r[i])))):.3f}")
print("total over forms with res > 0.01:", int(np.sum(per > 0.01)), "sum", per.sum())
