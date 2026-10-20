"""lm_alt.py (rm2-phase, restart session): rm2-solve's calibrated Levenberg-Marquardt (glm.py, copied unchanged) on
per-array targets that may be ALTERNATING (Lemma A).  The model only needs t_i = Q(z_i) and p_l = <Q, T_l> over all
binary quadratic forms Q, so nothing in it assumes T_l = h h^T.

usage: lm_alt.py m DIRS TARGETS ORDER restarts iters seed tlimit [tag]
  DIRS: comma list of directions (ints, bit j = coordinate j; multiplicity = repetition)
  TARGETS: comma list of target indices (t_small_alt.py order with ALT=1), or 'plain:<idx>,<idx>,...' = the plain
           program for those targets (DIRS ignored: the union of minimal decompositions)
  ORDER: lex | randN (random time order with seed N) | keep
Output: one line per restart, final cost on ALL forms and the per-array split; best x saved.
"""
import sys, time, os, warnings
import numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning)
from xbase import Forms
from glm import Model, lm, random_start
from tile_alt import target_matrices, decomposition


def build(m, dspec, tspec, order):
    mats, nonalt = target_matrices(m)
    if tspec.startswith("plain:"):
        tt = [int(x) for x in tspec[6:].split(",")]
        dirs = []
        for t in tt: dirs += decomposition(mats[t])
    else:
        tt = [int(x) for x in tspec.split(",")]
        dirs = [int(z) for z in dspec.split(",")]
    if order == "lex": dirs = sorted(dirs)
    elif order.startswith("rand"):
        r = np.random.default_rng(int(order[4:])); dirs = [dirs[i] for i in r.permutation(len(dirs))]
    fo = Forms(m)
    di = {z: i for i, z in enumerate(fo.dirs)}
    t = np.stack([fo.vals[:, di[z]] for z in dirs], axis=1).astype(float)          # forms x M
    p = np.stack([fo.parity(fo.tvec(mats[k])) for k in tt], axis=1).astype(float)    # forms x W
    keep = np.arange(1, fo.nf)                                                       # Q != 0
    return t[keep], p[keep], dirs, tt, nonalt


def main():
    a = sys.argv[1:]
    m = int(a[0]); dspec = a[1]; tspec = a[2]; order = a[3]
    restarts, iters, seed = int(a[4]), int(a[5]), int(a[6]); tlimit = float(a[7]); tag = a[8] if len(a) > 8 else "run"
    t, p, dirs, tt, nonalt = build(m, dspec, tspec, order)
    W, M = p.shape[1], t.shape[1]
    md = Model(t, p, W, M)
    print(f"[{tag}] m={m} W={W} M={M} targets {tt} ({sum(x >= nonalt for x in tt)} alternating) dirs(time) {dirs} "
          f"forms {t.shape[0]} unknowns {md.nvar}", flush=True)
    rng = np.random.default_rng(seed); best = np.inf
    for r in range(restarts):
        t0 = time.time()
        x0 = random_start(md, rng, 1.0, 0.3)
        x, c, hist = lm(md, x0, iters=iters, idx=None, tlimit=tlimit, verbose=0)
        call, pa = md.cost(x, None, per_array=True)
        print(f"[{tag}] RESTART {r} ALL {call:.6e} per-array {np.array2string(pa, precision=3)} its {len(hist)} "
              f"{time.time()-t0:.0f}s", flush=True)
        if call < best: best = call; np.save(f"best_{tag}.npy", x)
    print(f"[{tag}] BEST ALL {best:.6e}", flush=True)


if __name__ == "__main__":
    main()
