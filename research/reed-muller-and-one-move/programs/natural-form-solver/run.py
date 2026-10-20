"""run.py (rm2-solve): LM runs on one instance, several random restarts.

usage: run.py m TSPEC DIRS ORDER restarts iters nsub seed tlimit [sa sh] [tag]
  TSPEC: 'k:0,2,4,...' canonical Singer targets T_k; 'I:8' slice (all targets I)
  DIRS:  'all' | 'frames' | 'Iframes'
  ORDER: lex | singer | randN | keep
  nsub:  0 = all forms; else fit on a random subset of nsub forms (fixed per restart), reported on all forms
  tlimit: seconds per restart
Output: one line per restart: final cost on the subset and on ALL forms, per-array split; best x saved as .npy.
"""
import sys, time, warnings
import numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning)
from inst import canonical, slice_I
from glm import Model, lm, lm_weighted, random_start

a = sys.argv[1:]
m = int(a[0]); tspec = a[1]; dspec = a[2]; order = a[3]
restarts = int(a[4]); iters = int(a[5]); nsub = int(a[6]); seed = int(a[7]); tlimit = float(a[8])
sa = float(a[9]) if len(a) > 9 else 1.0; sh = float(a[10]) if len(a) > 10 else 0.3
tag = a[11] if len(a) > 11 else "run"
mode = a[12] if len(a) > 12 else "lm"
kind, rest = tspec.split(":")
if dspec.startswith("d:"): dspec = [int(v) for v in dspec[2:].split(",")]
if kind == "k":
    ks = [int(v) for v in rest.split(",")]
    inst, fd = canonical(m, ks, dirs=dspec, order=order, seed=seed)
elif kind == "I":
    inst, fd = slice_I(m, int(rest), dirs=dspec, order=order, seed=seed)
elif kind == "r":   # random survivor: 'r:W/s' = first random multiset of W canonical targets (seed s) with no tight form
    import random as _r
    Wr, sr = [int(v) for v in rest.split("/")]; rr = _r.Random(sr)
    while True:
        ks = sorted(rr.randrange(31) for _ in range(Wr))
        inst, fd = canonical(m, ks, dirs=dspec, order=order, seed=seed)
        et = inst.exact_tests()
        if et["span_ok"] and et["violations"] == 0 and et["tight"] == 0: break
    print(f"[{tag}] random survivor ks={ks}", flush=True)
print(f"[{tag}] {inst.label} W={inst.W} M={inst.M} forms={inst.Nf} exact={inst.exact_tests()}", flush=True)
print(f"[{tag}] dirs (time order) = {inst.dirs}", flush=True)
md = Model(inst.t, inst.p, inst.W, inst.M)
print(f"[{tag}] unknowns {md.nvar}", flush=True)
rng = np.random.default_rng(seed)
best = (np.inf, None)
for r in range(restarts):
    t0 = time.time()
    idx = None if nsub == 0 or nsub >= inst.Nf else np.sort(rng.choice(inst.Nf, nsub, replace=False))
    x0 = random_start(md, rng, sa, sh)
    if mode == "lm":
        x, c, hist = lm(md, x0, iters=iters, idx=idx, tlimit=tlimit, verbose=1, tag=f"[{tag} r{r}]")
    else:   # 'w': anti-split re-weighting, epochs of 30 LM iterations
        x, c, hist = lm_weighted(md, x0, epochs=max(1, iters // 30), per=30, idx=idx, tlimit=tlimit, tag=f"[{tag} r{r}]")
    call, pa = md.cost(x, None, per_array=True)
    A_, B_, H_ = md.unpack(x)
    print(f"[{tag}] norms |A| {np.abs(A_).max():.2e} |B| {np.abs(B_).max():.2e} |H| {np.abs(H_).max():.2e}", flush=True)
    print(f"[{tag}] RESTART {r} subset {c:.6e} ALL {call:.6e} per-array {np.array2string(pa, precision=3)} "
          f"its {len(hist)} {time.time() - t0:.0f}s", flush=True)
    if call < best[0]:
        best = (call, x); np.save(f"best_{tag}.npy", x)
print(f"[{tag}] BEST ALL {best[0]:.6e}", flush=True)
