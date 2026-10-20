"""homotopy.py (rm2-solve): deletion homotopy from a solved plain W-in-5W program to W in M_target.

usage: homotopy.py KS ORDERSEED NDEL DELSEED nsub mu0 mu1 steps per [x40.npy]
 1. plain instance (frames of the targets T_k, k in KS, random time order ORDERSEED); solve it from a random start
    (or load x40.npy);
 2. choose NDEL moves to delete (random among moves whose direction occurs more than once, DELSEED);
 3. penalty path: mu * (|in_i|^2 + |out_i|^2) on each deleted move (in_i = A[:,i], H[:i,i]; out_i = B[:,i], H[i,i+1:];
    the sum is minimised in the balanced gauge, so it measures the move's contribution, not its scale);
    mu = mu0 .. mu1 geometrically in `steps` steps, `per` LM iterations each; report the unpenalised residual;
 4. drop the deleted moves, polish the remaining (W in 5W - NDEL) program with plain LM on the subset; report on ALL forms.
Any directions are allowed: by Theorem B(a) a multi-target program with Singer targets is a no-phase Singer unit with
the same W and M, whatever its directions.
"""
import sys, time, warnings, random
import numpy as np
warnings.filterwarnings("ignore")
from inst import canonical, Instance
from glm import Model, lm, random_start

a = sys.argv[1:]
ks = [int(v) for v in a[0].split(",")]; oseed = int(a[1]); ndel = int(a[2]); dseed = int(a[3]); nsub = int(a[4])
mu0, mu1, steps, per = float(a[5]), float(a[6]), int(a[7]), int(a[8])
xfile = a[9] if len(a) > 9 else None
tag = f"hom k{len(ks)} o{oseed} d{ndel}s{dseed}"
import os
MM = int(os.environ.get("HM", "5"))
inst, fd = canonical(MM, ks, dirs="frames", order=f"rand{oseed}")
md = Model(inst.t, inst.p, inst.W, inst.M)
rng = np.random.default_rng(oseed * 1000 + dseed)
idx = np.sort(rng.choice(inst.Nf, nsub, replace=False)) if 0 < nsub < inst.Nf else None
print(f"[{tag}] {inst.label} M={inst.M} dirs={inst.dirs}", flush=True)
t0 = time.time()
if xfile:
    x = np.load(xfile)
else:
    for att in range(4):
        x, c, h = lm(md, random_start(md, rng, 1, 0.3), iters=400, idx=idx, verbose=0, tlimit=1500)
        print(f"[{tag}] plain solve attempt {att}: {c:.3e} its {len(h)} t {time.time() - t0:.0f}s", flush=True)
        if c < 1e-18: break
    np.save(f"x40_o{oseed}.npy", x)
M = inst.M
cnt = {}
for z in inst.dirs: cnt[z] = cnt.get(z, 0) + 1
dup = [i for i, z in enumerate(inst.dirs) if cnt[z] > 1]
r = random.Random(dseed)
# choose the deletion set among duplicates (never all copies of a direction) so that the resulting multiset passes
# span and deletion and has as few tight forms as possible (first one with none, else the best of the tries)
best_d = None
for tr in range(int(os.environ.get("HTRIES", "300"))):
    cand = dup[:]; r.shuffle(cand); left = dict(cnt); dd = []
    for i in cand:
        if len(dd) == ndel: break
        if left[inst.dirs[i]] > 1: dd.append(i); left[inst.dirs[i]] -= 1
    if len(dd) < ndel: continue
    kk = [i for i in range(M) if i not in dd]
    et = Instance(MM, inst.T, [inst.dirs[i] for i in kk]).exact_tests()
    if not et["span_ok"] or et["violations"] > 0: continue
    if best_d is None or et["tight"] < best_d[1]: best_d = (dd, et["tight"])
    if et["tight"] == 0: break
if best_d is None: print(f"[{tag}] no deletion set passes span+deletion", flush=True); sys.exit(0)
dele = best_d[0]
print(f"[{tag}] deletion set found after {tr + 1} tries, tight forms {best_d[1]}", flush=True)
dele = sorted(dele)
print(f"[{tag}] delete moves {dele} (directions {[inst.dirs[i] for i in dele]})", flush=True)
# penalty mask over packed variables
W = md.W; nA = md.nA; iu = md.iu
mask = np.zeros(md.nvar)
for i in dele:
    for c_ in range(W):
        mask[c_ * M + i] += 1          # A[c,i]
        mask[nA + c_ * M + i] += 1     # B[c,i]
hpos = {(int(p), int(q)): 2 * nA + k for k, (p, q) in enumerate(zip(*iu))}
for (p, q), k in hpos.items():
    if p in dele: mask[k] += 1
    if q in dele: mask[k] += 1
for s, mu in enumerate(np.geomspace(mu0, mu1, steps)):
    md.pen = mu * mask
    x, c, h = lm(md, x, iters=per, idx=idx, verbose=0, tlimit=3000)
    md.pen = None
    res = md.cost(x, idx); pn = float((mask * np.abs(x) ** 2).sum())
    print(f"[{tag}] mu {mu:.2e} residual {res:.4e} deleted-size {pn:.4e} t {time.time() - t0:.0f}s", flush=True)
# drop deleted moves
keep = [i for i in range(M) if i not in dele]
A, B, H = md.unpack(x)
inst2 = Instance(MM, inst.T, [inst.dirs[i] for i in keep], label=f"{inst.label} minus {dele}")
md2 = Model(inst2.t, inst2.p, inst2.W, inst2.M)
x2 = md2.pack(A[:, keep], B[:, keep], H[np.ix_(keep, keep)])
idx2 = np.sort(rng.choice(inst2.Nf, nsub, replace=False)) if 0 < nsub < inst2.Nf else None
print(f"[{tag}] after dropping: M={inst2.M} subset {md2.cost(x2, idx2):.4e} exact={inst2.exact_tests()}", flush=True)
x2, c2, h2 = lm(md2, x2, iters=300, idx=idx2, verbose=0, tlimit=2000)
call, pa = md2.cost(x2, None, per_array=True)
A2, B2, H2 = md2.unpack(x2)
print(f"[{tag}] FINAL W={inst2.W} M={inst2.M} subset {c2:.4e} ALL {call:.4e} per-array {np.array2string(pa, precision=3)} "
      f"norms {np.abs(A2).max():.1e} {np.abs(B2).max():.1e} {np.abs(H2).max():.1e} t {time.time() - t0:.0f}s", flush=True)
np.save(f"hom_o{oseed}_d{ndel}s{dseed}.npy", x2)
