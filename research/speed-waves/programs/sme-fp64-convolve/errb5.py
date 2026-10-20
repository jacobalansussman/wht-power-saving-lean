"""errb.py -- error of the blocked SME kernel (fc_blk tap-columns per block) vs numpy, against the exact sum."""
import sys
import numpy as np
import fc64
from check import exact_full, rng
for (n, m) in [(100000, 33), (100000, 257), (100000, 1025), (20000, 4097), (5000, 1000)]:
    a = rng.uniform(-1, 1, n); v = rng.uniform(-1, 1, m)
    for mode in ("same", "valid"):
        o0, nout = fc64._range(n, m, mode)
        ts = sorted(set(list(range(min(16, nout))) + list(range(max(0, nout - 16), nout)) + rng.integers(0, nout, 300).tolist()))
        ex, sc = exact_full(a, v[::-1], [o0 + t for t in ts])
        ynp = np.convolve(a, v, mode)
        enp = np.abs(ynp[ts] - ex) / sc
        line = f"n={n:6d} M={m:5d} {mode:5s} numpy max {enp.max():.2e} med {np.median(enp):.2e} |"
        for blk in (0, 128, 64, 32):
            fc64.set_global("fc_blk5", blk if blk else 100000)
            y = fc64.convolve(a, v, mode, path=5)
            d = np.max(np.abs(y - ynp)) / np.max(np.abs(ynp))
            e = np.abs(y[ts] - ex) / sc
            line += f" blk{blk}: max {e.max():.2e} med {np.median(e):.2e} (vs np {d:.1e}) |"
        print(line); sys.stdout.flush()
fc64.set_global("fc_blk5", 64)
