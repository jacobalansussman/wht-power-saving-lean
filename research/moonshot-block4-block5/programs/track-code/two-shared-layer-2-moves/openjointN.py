"""openjointN.py -- w5-twoshared: joint test of the Lemma-D-open cores, NUMERICAL first (the exact systems have ~100
unknowns): for pairs (top-0 set, top-1 set) of Lemma-D injections, LM over all program coefficients of the summed
closure residual of BOTH outputs (readouts by least squares), several starts; a pair with residual < 1e-12 is saved
for the exact decision (Singular) and verify.py; the others are reported with their floors.  Floors are local minima
(no proof); the control is ctl_verify (plain 2 x 4 reaches 5e-32 from random starts).
Usage: openjointN.py seed core npairs starts cap"""
import sys, os, itertools, pickle, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
seed, cidx, npairs, starts, cap = [int(a) for a in sys.argv[1:6]]
sys.argv = ['opencase.py', str(seed), str(cidx), '1', str(cap)]
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'opencase.py')).read().split("# (i) each output alone")[0]
exec(compile(src, 'opencase_build', 'exec'))
from exact import Closure
from scipy.optimize import least_squares
pairs = list(itertools.product(range(len(inj_by_top[0])), range(len(inj_by_top[1]))))
random.Random(5).shuffle(pairs)
nr = np.random.default_rng(7)
floors = []
for (a, b) in pairs[:npairs]:
    mvs = moves_of(R0, R1, None, None, list(inj_by_top[0][a][1]) + list(inj_by_top[1][b][1]))
    P = build(2, mvs)
    Cs = [Closure(P, l, cap) for l in range(2)]
    if not all(C.ok for C in Cs):
        print(f"pair ({a},{b}): closure cap", flush=True); continue
    k = P.np
    def fun(x):
        par = x[:k] + 1j * x[k:]
        v = np.concatenate([C.resid_vec(par) for C in Cs])
        return np.concatenate([v.real, v.imag])
    best = np.inf
    for s in range(starts):
        r = least_squares(fun, nr.standard_normal(2 * k), method='trf', max_nfev=1500)
        best = min(best, 2 * r.cost)
    floors.append(best)
    print(f"pair ({a},{b}) types {inj_by_top[0][a][0]} / {inj_by_top[1][b][0]}: closures {[(C.nr, C.nc) for C in Cs]}, floor {best:.4g}", flush=True)
    if best < 1e-12:
        fn = f'open_sing/NUMHIT_s{seed}_c{cidx}_{a}_{b}.pkl'
        pickle.dump({'moves': mvs}, open(fn, 'wb'))
        print("   NUMERICAL HIT saved", fn, flush=True)
print("floors sorted", np.round(sorted(floors), 4).tolist(), flush=True)
