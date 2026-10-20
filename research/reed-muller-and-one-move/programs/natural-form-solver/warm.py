"""warm.py KS ORDER SEED ITERS NSUB TLIMIT [eps]: 8 in 31 (one move per direction) started from the LARGEST plain
packing: the arrays of a maximum family of pairwise disjoint frames (packing2.py) get the exact XOR-chain program
(a = b = 1 on their frame, H = 2 between moves of the same frame; F = sum over chains (-2)^(k-1) t^C = XOR), the
other arrays and all other entries start at eps x random. Then LM. Reports which arrays end served."""
import sys, time, warnings; warnings.filterwarnings("ignore")
import numpy as np
from inst import canonical
from glm import Model, lm
from packing2 import packing
ks = [int(v) for v in sys.argv[1].split(",")]; order = sys.argv[2]; seed = int(sys.argv[3])
iters, nsub, tlimit = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]); eps = float(sys.argv[7]) if len(sys.argv) > 7 else 0.1
inst, fd = canonical(5, ks, dirs="all", order=order)
md = Model(inst.t, inst.p, inst.W, inst.M); W, M = inst.W, inst.M
pn, fam = packing(ks)
rng = np.random.default_rng(seed)
def cn(*s): return (rng.standard_normal(s) + 1j * rng.standard_normal(s)) / np.sqrt(2)
A = eps * cn(W, M); B = eps * cn(W, M); H = np.triu(eps * cn(M, M), 1)
pos = {z: i for i, z in enumerate(inst.dirs)}
packed = [l for l, F in fam]
for l, F in fam:
    idx = sorted(pos[z] for z in F)
    A[:, idx] = 0; B[:, idx] = 0; A[l, idx] = 1; B[l, idx] = 1
    for a_ in idx:
        H[a_, :] = 0; H[:, a_] = 0
    for p_, a_ in enumerate(idx):
        for b_ in idx[p_ + 1:]: H[a_, b_] = 2
x0 = md.pack(A, B, H)
c0, pa0 = md.cost(x0, per_array=True)
print(f"[warm] ks={ks} order={order} packing {pn} arrays {packed}; start ALL {c0:.4e} per-array {np.array2string(pa0, precision=2)}", flush=True)
idx = np.sort(rng.choice(inst.Nf, nsub, replace=False)) if 0 < nsub < inst.Nf else None
t0 = time.time()
x, c, h = lm(md, x0, iters=iters, idx=idx, verbose=0, tlimit=tlimit)
call, pa = md.cost(x, per_array=True)
print(f"[warm] END ALL {call:.4e} per-array {np.array2string(pa, precision=3)} served {[l for l in range(W) if pa[l] < 1e-10]} its {len(h)} {time.time() - t0:.0f}s", flush=True)
np.save(f"warm_{order}_{seed}.npy", x)
