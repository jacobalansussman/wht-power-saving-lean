# check_big.py -- long segments (ObsPy PPSD-like nperseg 2^13..2^15) against scipy.signal.welch, every bin
import sys, os, numpy as np, scipy.signal
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import as0py
rng = np.random.default_rng(5); nf = 0
for N, eng in ((8192, "smed"), (16384, "smed"), (32768, "smef"), (8192, "kfr"), (32768, "kfr"), (32768, "vdsp")):
    for L in (N * 9 + 7, N * 12):
        x = rng.standard_normal(L) + 0.5
        r = scipy.signal.welch(x, fs=100.0, window="hann", nperseg=N)[1]
        a = as0py.welch(x, fs=100.0, window="hann", nperseg=N, engine=eng)[1]
        e = np.max(np.abs(a - r) / np.abs(r)); nf += e >= 1e-13
        print(f"CHECK N={N} {eng} L={L} nseg={(L - N) // (N // 2) + 1} max bin rel {e:.2e} {'PASS' if e < 1e-13 else 'FAIL'}")
print("FAILS", nf)
