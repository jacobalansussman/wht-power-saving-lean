"""tri_direct.py CASES_FILE NRESTART ITERS NSUB: random-start LM on the multisets 'plain frames minus one move' that
pass the exact tests (lines 'PASS ks=(..) frames=(..) delete dir z ...' of scan_frames.py). Each case: a random
time order (seed = case number), NRESTART restarts; residual on ALL forms; per-array split."""
import sys, re, time, warnings; warnings.filterwarnings("ignore")
import numpy as np
from inst import canonical
from glm import Model, lm, random_start
from packing2 import all_frames
cases = [l for l in open(sys.argv[1]) if l.startswith("PASS")]
nres, iters, nsub = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
for ci, line in enumerate(cases):
    m = re.search(r"ks=\(([0-9, ]+)\) frames=\(([0-9, ]+)\) delete dir (\d+)", line)
    ks = [int(v) for v in m.group(1).split(",")]; ch = [int(v) for v in m.group(2).split(",")]; z = int(m.group(3))
    d = [w for k, c in zip(ks, ch) for w in sorted(sorted(F) for F in all_frames(k))[c]]
    d.remove(z)
    inst, fd = canonical(int(__import__("os").environ.get("HM", "5")), ks, dirs=d, order=f"rand{100 + ci}")
    md = Model(inst.t, inst.p, inst.W, inst.M); rng = np.random.default_rng(ci)
    best = np.inf
    for r in range(nres):
        idx = np.sort(rng.choice(inst.Nf, nsub, replace=False)) if 0 < nsub < inst.Nf else None
        t0 = time.time()
        x, c, h = lm(md, random_start(md, rng, 1, 0.3), iters=iters, idx=idx, verbose=0, tlimit=900)
        call, pa = md.cost(x, per_array=True); best = min(best, call)
        print(f"case {ci} ks={ks} frames={ch} minus {z}: W={inst.W} M={inst.M} exact={inst.exact_tests()} restart {r}: "
              f"ALL {call:.4e} per-array {np.array2string(pa, precision=2)} its {len(h)} {time.time() - t0:.0f}s", flush=True)
    print(f"CASE {ci} BEST {best:.4e}", flush=True)
