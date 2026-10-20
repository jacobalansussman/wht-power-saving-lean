# au-f: is as-0's Python baseline (scipy.signal.welch) the strongest a practitioner would write?
# A 6-line numpy/scipy.fft Welch (batched rfft, workers) in fp64 and fp32, against scipy.signal.welch and as-0's SME
# engine (ctypes), 64 channels x 2^18 float64, 50 % Hann. Programs take turns; 9 repeats of >= 0.2 s; ms per call.
import os, sys, time, numpy as np, scipy.signal, scipy.fft
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import as0py
from numpy.lib.stride_tricks import sliding_window_view as swv
rng = np.random.default_rng(7); fs = 1000.0
x = rng.standard_normal((64, 1 << 18)) + 0.3
def fast(x, N, workers, dt):
    step = N // 2; w = scipy.signal.get_window('hann', N).astype(dt)
    seg = swv(x.astype(dt, copy=False), N, axis=-1)[..., ::step, :]
    seg = (seg - seg.mean(-1, keepdims=True)) * w
    X = scipy.fft.rfft(seg, axis=-1, workers=workers)
    P = (X.real**2 + X.imag**2).mean(axis=-2) * (1.0 / (fs * (w.astype(np.float64)**2).sum()))
    P[..., 1:-1] *= 2
    return P
for N in [int(a) for a in sys.argv[1:]] or [512, 4096]:
    x32 = x.astype(np.float32)
    progs = [("scipy.signal.welch 1 thread", lambda: scipy.signal.welch(x, fs=fs, nperseg=N)[1]),
             ("numpy 6-liner fp64 rfft workers=1", lambda: fast(x, N, 1, np.float64)),
             ("numpy 6-liner fp64 rfft workers=4", lambda: fast(x, N, 4, np.float64)),
             ("numpy 6-liner fp32 rfft workers=4", lambda: fast(x32, N, 4, np.float32)),
             ("as-0 SME engine (ctypes, 1 thread)", lambda: as0py.welch(x, fs=fs, nperseg=N, engine="sme")[1])]
    ref = scipy.signal.welch(x, fs=fs, nperseg=N)[1]
    for name, f in progs:
        e = np.abs(f() - ref).max() / np.abs(ref).max()
        print(f"CHECK N={N} {name:38s} max rel (vs max bin) {e:.2e}")
    ks = []
    for name, f in progs:
        k = 1
        while True:
            t0 = time.perf_counter(); [f() for _ in range(k)]
            if time.perf_counter() - t0 >= 0.2: break
            k *= 2
        ks.append(k)
    T = [[] for _ in progs]
    for rep in range(9):
        order = list(range(len(progs))); order = order[rep % len(order):] + order[:rep % len(order)]
        for i in order:
            t0 = time.perf_counter(); [progs[i][1]() for _ in range(ks[i])]; T[i].append((time.perf_counter() - t0) / ks[i] * 1e3)
    for i, (name, f) in enumerate(progs):
        t = sorted(T[i]); print(f"RESULT N={N} {name:38s} best {t[0]:8.2f} median {t[4]:8.2f} ms  reps {' '.join(f'{v:.1f}' for v in T[i])}", flush=True)
