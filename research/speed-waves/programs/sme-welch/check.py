# check.py -- every engine of libas0 against scipy.signal.welch, max relative error over EVERY bin
# (max |a - r| / max |r| and also max over bins of |a - r| / |r|), odd and even segment counts, two windows,
# two overlaps, plus deliberately broken copies (AS0_BROKEN=1 skips the fold) that must FAIL.
import sys, os, numpy as np, scipy, scipy.signal
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import as0py
broken = os.environ.get("AS0_BROKEN", "0") != "0"
rng = np.random.default_rng(11)
worst = 0.0; nfail = 0
for N in (512, 1024, 2048, 4096):
    for window, nov in (("hann", N // 2), ("hamming", 0), ("hann", N // 4)):
        for L in (N * 37 + 5, N * 40, (1 << 18) + 3):   # segment counts odd and even
            x = rng.standard_normal(L) * 2.0 + 0.7 + 0.3 * np.sin(0.05 * np.arange(L))   # DC offset: detrend matters
            f, r = scipy.signal.welch(x, fs=250.0, window=window, nperseg=N, noverlap=nov)
            nseg = (L - N) // (N - nov) + 1
            for eng in ("sme", "smef", "smed", "vdsp", "pfd", "kfr", "smew", "fftw"):
                if eng in ("smef", "smew") and N not in (512, 4096): continue
                if eng == "smed" and N not in (1024, 2048): continue
                _, a = as0py.welch(x, fs=250.0, window=window, nperseg=N, noverlap=nov, engine=eng)
                e_glob = np.max(np.abs(a - r)) / np.max(np.abs(r)); e_bin = np.max(np.abs(a - r) / np.abs(r))
                ok = e_bin < 1e-13
                if not broken: worst = max(worst, e_bin)
                nfail += (not ok)
                print(f"CHECK N={N:5d} {window:7s} nov={nov:5d} L={L:7d} nseg={nseg:5d} ({'odd ' if nseg % 2 else 'even'}) {eng:5s} "
                      f"max|d|/max|r| {e_glob:.2e}  max bin rel {e_bin:.2e}  {'PASS' if ok else 'FAIL'}")
# multichannel (2-D) path
x = rng.standard_normal((5, 100003)); f, r = scipy.signal.welch(x, fs=1000.0, window="hamming", nperseg=1024, noverlap=0)
for eng in ("sme", "smed", "vdsp", "pfd", "kfr"):
    _, a = as0py.welch(x, fs=1000.0, window="hamming", nperseg=1024, noverlap=0, engine=eng)
    e = np.max(np.abs(a - r) / np.abs(r)); nfail += e >= 1e-13
    print(f"CHECK 2-D 5 x 100003 hamming nov=0 N=1024 {eng:5s} max bin rel {e:.2e} {'PASS' if e < 1e-13 else 'FAIL'}")
print(f"numpy {np.__version__} scipy {scipy.__version__}; broken={broken}; FAILS {nfail}; worst bin rel (unbroken) {worst:.2e}")
