# cal1.py -- calibration on a2-nag's exported frame-sum relaxations: existence (plain cost) and refusal (2 in 5, border)
import sys, json, time, glob, numpy as np
from certnum import *
out = open(sys.argv[2], 'w') if len(sys.argv) > 2 else sys.stdout
def log(*a):
    print(*a, file=out, flush=True)
files = sorted(glob.glob(sys.argv[1]))
rng = np.random.default_rng(1)
for fn in files:
    t0 = time.time(); F = PolySys.from_a2nag_json(fn); info = F.info or {}
    na = info.get('nalpha_free', 0)
    best = None
    for st in range(12):
        x0 = (rng.normal(size=F.n) + 1j * rng.normal(size=F.n))
        x, f = lm_solve(F, x0, iters=300)
        if best is None or f < best[1]: best = (x, f)
        if f < 1e-26: break
    x, f = best
    log('=== %s  n=%d N=%d  nalpha=%d  LM best |F|^2=%.3e  max|x|=%.3g' % (fn.split('/')[-1], F.n, F.N, na, f, np.abs(x).max()))
    if f < 1e-20:
        rep = certify_point(F, x, prefer_fix=list(range(na)))
        for k in ['rank', 'sv_gap', 'square', 'fixed', 'krawczyk', 'rho', 'leftover_rows', 'leftover_max_abs_upper',
                  'exact_candidate_maxden', 'exact_nonzero_rows', 'verdict', 'exact_point']:
            if k in rep: log('   %s: %s' % (k, rep[k]))
    else:
        r = F.eval_np(x); y = np.conj(r)
        ok, yq, const = dual_from_point(F, x)
        log('   dual certificate (degree 0) from the near miss: %s%s' % (ok, ('  constant=%s, support %d rows' % (const, sum(1 for c in yq if c != GZERO))) if ok else ''))
        rho, how = exclusion_radius(F, x, y)
        log('   certified exclusion box around the near miss: half-width %s (%s)' % (rho, how))
        for d in [0, 1, 2]:
            t1 = time.time(); found, msg = nulla(F, d)
            log('   NulLA degree %d (exact over Q): %s  [%s, %.1fs]' % (d, found, msg, time.time() - t1))
            if found: break
    log('   time %.1fs' % (time.time() - t0))
