"""served.py SPEC ORDER XFILE: restrict a solved run to its exactly served arrays; for each move, the residual of
that sub-program when the move is dropped (no re-fit); the essential moves of each served array."""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np
from inst import canonical, Instance
from glm import Model
from packing import frames
spec, order, xf = sys.argv[1:4]
ks = [int(v) for v in spec.split(",")]
inst, fd = canonical(5, ks, dirs="all", order=order)
md = Model(inst.t, inst.p, inst.W, inst.M)
x = np.load(xf)
c, pa = md.cost(x, per_array=True)
served = [l for l in range(inst.W) if pa[l] < 1e-10]
print("ALL", c, "served", served, "ks", [ks[l] for l in served])
A, B, H = md.unpack(x)
def sub(arrs, keep):
    it = Instance(5, [inst.T[l] for l in arrs], [inst.dirs[i] for i in keep])
    m2 = Model(it.t, it.p, it.W, it.M)
    x2 = m2.pack(A[np.ix_(arrs, keep)], B[np.ix_(arrs, keep)], H[np.ix_(keep, keep)])
    return m2.cost(x2, per_array=True)
M = inst.M
c0, p0 = sub(served, list(range(M)))
print("sub-program residual", c0)
ess = {l: [] for l in served}
for i in range(M):
    keep = [j for j in range(M) if j != i]
    ci, pi = sub(served, keep)
    hit = [served[q] for q in range(len(served)) if pi[q] > 1e-8]
    for l in hit: ess[l].append(i)
    print(f"move {i:2d} dir {inst.dirs[i]:2d}: drop -> residual {ci:.3e}, breaks arrays {hit}")
fr = frames(ks)
for l in served:
    dset = sorted(inst.dirs[i] for i in ess[l])
    print(f"array {l} (k={ks[l]}): essential moves {len(ess[l])}: dirs {dset}; frame {sorted(fr[l])}")
allm = sorted(set(i for l in served for i in ess[l]))
print("union of essential moves:", len(allm), "for", len(served), "arrays (plain would need", 5 * len(served), ")")
