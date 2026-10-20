"""lm_batch_alt.py (rm2-phase, restart session): Levenberg-Marquardt (rm2-solve's glm.py) on every class that the
tiling test lets through, read from tile_alt.py's output lines "class i: dirs (...) targets (...) ...: passes ...".
Per class: several random time orders x restarts; plus the calibration pair first (plain 4 in 18 with two alternating
targets must be found; a zero-block-excluded 3 in 11 must be refused).

usage: python3 lm_batch_alt.py m tile_output.txt [maxclasses] [orders] [restarts] [iters] [tlimit]
"""
import sys, re, time
import numpy as np
from glm import Model, lm, random_start
from lm_alt import build


def run(m, dspec, tspec, order, restarts, iters, tlimit, seed, tag):
    t, p, dirs, tt, nonalt = build(m, dspec, tspec, order)
    W, M = p.shape[1], t.shape[1]
    md = Model(t, p, W, M); rng = np.random.default_rng(seed); best = (np.inf, None)
    for r in range(restarts):
        x0 = random_start(md, rng, 1.0, 0.3)
        x, c, hist = lm(md, x0, iters=iters, idx=None, tlimit=tlimit, verbose=0)
        call, pa = md.cost(x, None, per_array=True)
        if call < best[0]: best = (call, pa)
    print(f"  [{tag}] order {order}: best ALL {best[0]:.4e} per-array {np.array2string(best[1], precision=2)}",
          flush=True)
    return best[0]


def main():
    a = sys.argv[1:]
    m = int(a[0]); path = a[1]
    maxcl = int(a[2]) if len(a) > 2 else 40; norders = int(a[3]) if len(a) > 3 else 4
    restarts = int(a[4]) if len(a) > 4 else 2; iters = int(a[5]) if len(a) > 5 else 300
    tlimit = float(a[6]) if len(a) > 6 else 60
    t0 = time.time()
    print("calibration:", flush=True)
    run(4, "-", "plain:428,440,0,7", "rand1", 2, iters, tlimit, 1, "pos plain 4 in 18, 2 alternating")
    run(4, "1,1,2,2,4,5,6,8,11,12,15", "254,355,442", "rand1", 2, iters, tlimit, 3, "neg 3 in 11 (Z-excluded)")
    pat = re.compile(r"class (\d+): dirs \(([\d, ]+)\) targets \(([\d, ]+?),?\) .*: passes")
    cls = []
    for L in open(path):
        mm = pat.search(L)
        if mm: cls.append((int(mm.group(1)), mm.group(2).replace(" ", ""), mm.group(3).replace(" ", "")))
    print(f"{len(cls)} passing classes in {path}; running {min(maxcl, len(cls))}", flush=True)
    floors = []
    for ci, d, tt in cls[:maxcl]:
        b = min(run(m, d, tt, f"rand{o}", restarts, iters, tlimit, 100 * ci + o, f"class {ci} targets {tt}")
                for o in range(1, norders + 1))
        floors.append((b, ci, tt))
    floors.sort()
    print("lowest floors:", [(f"{b:.3e}", ci, tt) for b, ci, tt in floors[:10]])
    print(f"done in {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
