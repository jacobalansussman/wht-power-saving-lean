"""check_auto.py -- the dispatching fc64.convolve / correlate (no path given) against numpy, every route."""
import numpy as np
import fc64
rng = np.random.default_rng(5)
res = []
for (n, m) in [(256, 3), (256, 129), (1024, 33), (5000, 9), (100000, 3), (20000, 257), (30000, 3000), (12000, 4097), (7, 300)]:
    for cplx in (False, True):
        a = rng.uniform(-1, 1, n) + (1j * rng.uniform(-1, 1, n) if cplx else 0)
        v = rng.uniform(-1, 1, m) + (1j * rng.uniform(-1, 1, m) if cplx else 0)
        for mode in ("full", "same", "valid"):
            o0, nout = fc64._range(max(n, m), min(n, m), mode)
            route = fc64._choose(nout, min(n, m), cplx)
            r = np.convolve(a, v, mode); g = fc64.convolve(a, v, mode)
            d = np.max(np.abs(g - r)) / np.max(np.abs(r))
            ok = g.shape == r.shape and g.dtype == r.dtype and d < 1e-12
            res.append(ok)
            print(f"{'PASS' if ok else 'FAIL'} n={n} m={m} {'c128' if cplx else 'f64'} {mode} route={route} rel diff {d:.1e}")
print(f"SUMMARY {sum(res)} of {len(res)} PASS")
