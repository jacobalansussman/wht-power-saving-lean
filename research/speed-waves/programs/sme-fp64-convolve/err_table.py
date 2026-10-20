"""err_table.py -- error of every route against the correctly rounded direct sum (av-1, wave 15).
Error of an output = |got - exact| / sum_k |g_k x_{t+k}| ("scale"); we report max and median over sampled outputs
(first/last 16 outputs plus 300 random ones), for 'same' mode convolutions of uniform [-1,1) data.
Also the max error relative to max|y| (the norm-wise view FFT users quote)."""
import sys

import numpy as np
import scipy.signal as ss

import fc64
from check import exact_full, rng

rows = []
cases = [(100000, 3), (100000, 33), (100000, 257), (100000, 1025), (20000, 4097)]
for cplx in (False, True):
    for (n, m) in cases:
        if cplx:
            a = rng.uniform(-1, 1, n) + 1j * rng.uniform(-1, 1, n)
            v = rng.uniform(-1, 1, m) + 1j * rng.uniform(-1, 1, m)
        else:
            a = rng.uniform(-1, 1, n)
            v = rng.uniform(-1, 1, m)
        mode = "same"
        o0 = fc64._range(n, m, mode)[0]
        ts = sorted(set(list(range(16)) + list(range(n - 16, n)) + rng.integers(0, n, 300).tolist()))
        ex, sc = exact_full(a, v[::-1], [o0 + t for t in ts])
        R = {"fc_sme": fc64.convolve(a, v, mode, path=1), "numpy": np.convolve(a, v, mode),
             "sp_oa": ss.oaconvolve(a, v, mode), "sp_fft": ss.fftconvolve(a, v, mode)}
        if not cplx:
            R["fc_neon"] = fc64.convolve(a, v, mode, path=2)
            R["vdsp"] = fc64.vdsp_convolve(a, v, mode)
        ymax = np.max(np.abs(ex))
        for k, y in R.items():
            e = np.abs(y[ts] - ex) / sc
            rows.append((("c128" if cplx else "f64"), n, m, k, np.max(e), np.median(e), np.max(np.abs(y[ts] - ex)) / ymax))
            print(f"{rows[-1][0]:4s} n={n:6d} M={m:5d} {k:8s} max {rows[-1][4]:.2e} median {rows[-1][5]:.2e}   "
                  f"max|err|/max|y| {rows[-1][6]:.2e}")
        sys.stdout.flush()
