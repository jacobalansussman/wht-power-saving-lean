"""warm2.py KS ORDER SEED ITERS NSUB TLIMIT [eps]: 8 in 31 (one move per direction) started with EVERY array on an
XOR chain over one of its frames, frames chosen to minimise the number of moves used by two or more arrays
(overlaps are forced: 40 frame slots on 31 moves). a = b = 1 on each array's frame, H = 2 between moves of a common
frame; everything else eps x random. The cross terms on shared moves are the residual to remove. Then LM."""
import sys, time, itertools, warnings; warnings.filterwarnings("ignore")
import numpy as np
from inst import canonical
from glm import Model, lm
from packing2 import all_frames
ks = [int(v) for v in sys.argv[1].split(",")]; order = sys.argv[2]; seed = int(sys.argv[3])
iters, nsub, tlimit = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]); eps = float(sys.argv[7]) if len(sys.argv) > 7 else 0.05
inst, fd = canonical(5, ks, dirs="all", order=order)
md = Model(inst.t, inst.p, inst.W, inst.M); W, M = inst.W, inst.M
frs = [sorted(sorted(F) for F in all_frames(k)) for k in ks]
rng = np.random.default_rng(seed)
best = None
for ch in (tuple(rng.integers(0, 6, W)) for _ in range(20000)):   # random sample of the 6^W frame choices
    cnt = np.zeros(32, dtype=int)
    for l in range(W):
        for z in frs[l][ch[l]]: cnt[z] += 1
    over = int((cnt[cnt > 1] - 1).sum()); shared = int((cnt > 1).sum())
    key = (over, shared, rng.random())
    if best is None or key < best[0]: best = (key, ch)
ch = best[1]
print(f"[warm2] ks={ks} order={order} frames {ch}: extra incidences {best[0][0]}, shared moves {best[0][1]}", flush=True)
def cn(*s): return (rng.standard_normal(s) + 1j * rng.standard_normal(s)) / np.sqrt(2)
A = eps * cn(W, M); B = eps * cn(W, M); H = np.triu(eps * cn(M, M), 1)
pos = {z: i for i, z in enumerate(inst.dirs)}
for l in range(W):
    idx = sorted(pos[z] for z in frs[l][ch[l]])
    A[l, idx] = 1; B[l, idx] = 1
    for p_, a_ in enumerate(idx):
        for b_ in idx[p_ + 1:]: H[a_, b_] = 2
x0 = md.pack(A, B, H)
c0, pa0 = md.cost(x0, per_array=True)
print(f"[warm2] start ALL {c0:.4e} per-array {np.array2string(pa0, precision=1)}", flush=True)
idx = np.sort(rng.choice(inst.Nf, nsub, replace=False)) if 0 < nsub < inst.Nf else None
t0 = time.time()
x, c, h = lm(md, x0, iters=iters, idx=idx, verbose=0, tlimit=tlimit)
call, pa = md.cost(x, per_array=True)
print(f"[warm2] END ALL {call:.4e} per-array {np.array2string(pa, precision=3)} served {[l for l in range(W) if pa[l] < 1e-10]} its {len(h)} {time.time() - t0:.0f}s", flush=True)
np.save(f"warm2_{order}_{seed}.npy", x)
