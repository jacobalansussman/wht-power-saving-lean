# calS.py -- Route S calibration: slice so that the root is IRRATIONAL (fix beta-type unknowns, solve for the alphas),
# so the simplest-rational closure cannot work; Route S must still certify existence on the solvable relaxations and
# must refuse (or fail to close) on the infeasible ones.
import sys, glob, json, numpy as np
from certnum import *
rng = np.random.default_rng(11)
for fn in sorted(glob.glob(sys.argv[1])):
    F = PolySys.from_a2nag_json(fn); na = (F.info or {}).get('nalpha_free', 0)
    best = None
    for st in range(12):
        x, f = lm_solve(F, rng.normal(size=F.n) + 1j * rng.normal(size=F.n), iters=300)
        if best is None or f < best[1]: best = (x, f)
        if f < 1e-26: break
    x, f = best
    print('=== %s n=%d N=%d LM %.1e' % (fn.split('/')[-1], F.n, F.N, f), flush=True)
    if f > 1e-20: print('   no numerical solution; skipped'); continue
    J = F.jac_np(x); r, _ = numrank(J)
    # choose the slice among the NON-alpha unknowns that are nonzero at x (so the alphas stay determined)
    rows = pivoted_qr_cols(J.T, r); Jc = J[rows, :].copy(); Jc[:, na:] *= 1e-3          # prefer alphas as determined
    cols = pivoted_qr_cols(Jc, r); fixv = [v for v in range(F.n) if v not in cols]
    mode = sys.argv[2] if len(sys.argv) > 2 else 'coord'
    rep = certify_point_S(F, x, fix_vars=fixv if mode == 'coord' else None, slice_mode=mode, rng=np.random.default_rng(5), degrees=tuple(tuple(int(c) for c in t.split(',')) for t in (sys.argv[3] if len(sys.argv) > 3 else '0,0;0,1;1,1').split(';')))
    for k in ['rank', 'slices', 'square', 'krawczyk', 'leftover', 'simplest_rational_closes', 'root_sample', 'route_S_closed', 'verdict']:
        if k in rep: print('   %s: %s' % (k, rep[k]), flush=True)
    print('   sliced unknowns:', [F.names[v] for v in fixv], ' alphas:', na, flush=True)
