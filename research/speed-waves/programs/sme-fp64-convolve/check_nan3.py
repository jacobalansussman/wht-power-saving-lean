"""check_nan.py -- non-finite inputs: the NaN / inf pattern of the output must equal numpy's direct sum exactly
(an FFT route would turn the whole output into NaN). Run with fc_nancheck = 1 (default) and = 0 (shows the problem)."""
import sys
import numpy as np
import fc64
rng = np.random.default_rng(11)
res = []
for nc in (3, 0):
    fc64.set_global("fc_nancheck", nc) if False else __import__("ctypes").c_int.in_dll(fc64._glib(), "fc_nancheck").__setattr__("value", nc)
    for (n, m) in ((5000, 33), (100000, 257), (20000, 4)):
        for cplx in (False, True):
            for bad in (np.nan, np.inf, -np.inf):
                a = rng.uniform(-1, 1, n) + (1j * rng.uniform(-1, 1, n) if cplx else 0)
                v = rng.uniform(-1, 1, m) + (1j * rng.uniform(-1, 1, m) if cplx else 0)
                pos = rng.integers(0, n, 3)
                a[pos] = bad
                for mode in ("same", "full", "valid"):
                    r = np.convolve(a, v, mode)
                    g = fc64.convolve(a, v, mode, path=1)
                    with np.errstate(invalid="ignore"):
                        same_pattern = np.array_equal(np.isnan(r), np.isnan(g)) and np.array_equal(np.isinf(r), np.isinf(g))
                        fin = np.isfinite(r)
                        d = np.max(np.abs(g[fin] - r[fin])) / np.max(np.abs(r[fin])) if fin.any() else 0.0
                    ok = same_pattern and d < 1e-12
                    if nc:
                        res.append(ok)
                    print(f"nancheck={nc} {'PASS' if ok else 'FAIL'} n={n} m={m} {'c128' if cplx else 'f64'} bad={bad} {mode}: "
                          f"nonfinite numpy {np.sum(~np.isfinite(r))} fc64 {np.sum(~np.isfinite(g))} finite rel diff {d:.1e}")
print(f"SUMMARY (nancheck=3) {sum(res)} of {len(res)} PASS")
