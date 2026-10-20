# aq-cross hybrid H1: scipy.signal.welch (float64) with ap-0's SME fp64 complex FFT and am-2/ap-2's
# "two real frames per complex FFT" packing. Whole call against scipy.signal.welch (scipy 1.13.1, pocketfft),
# plus a control with the SAME numpy glue but scipy.fft.fft (pocketfft) instead of the SME call.
# usage: python welch_sme.py NPERSEG LOG2LEN [check|time]
import ctypes, sys, time, os
import numpy as np, scipy, scipy.signal, scipy.fft
from numpy.lib.stride_tricks import sliding_window_view
AP0 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sme-fp64-batched-fft")
lib = ctypes.CDLL(os.path.join(AP0, "libap0v2.dylib"))
lib.ap0_new_v.restype = ctypes.c_void_p; lib.ap0_new_v.argtypes = [ctypes.c_long]
lib.ap0_fft_v.restype = None; lib.ap0_fft_v.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long]

N = int(sys.argv[1]); L = 1 << int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else "time"
fs = 1.0
win = scipy.signal.get_window('hann', N)
scale = 1.0 / (fs * (win * win).sum())
step = N // 2
nseg = (L - N) // step + 1
npair = nseg // 2
h = lib.ap0_new_v(N)
z = np.empty((npair, N), np.complex128)
Z = np.empty((npair, N), np.complex128)
neg = (-np.arange(N)) % N

def glue_pack(x):
    fr = sliding_window_view(x, N)[::step][:nseg]
    m = fr.mean(axis=1)
    zr = z.real; zi = z.imag
    np.subtract(fr[0:2 * npair:2], m[0:2 * npair:2, None], out=zr); zr *= win
    np.subtract(fr[1:2 * npair:2], m[1:2 * npair:2, None], out=zi); zi *= win
    return fr, m

def finish(P, fr, m):
    S = 0.5 * (P + P[neg])
    if nseg % 2:   # odd last segment: one real FFT
        last = (fr[nseg - 1] - m[nseg - 1]) * win
        X = np.fft.fft(last); S = S + (X.real ** 2 + X.imag ** 2)
    psd = S[:N // 2 + 1] * (scale / nseg)
    psd[1:-1] *= 2.0   # N even: DC and Nyquist not doubled
    return psd

def welch_sme(x):
    fr, m = glue_pack(x)
    lib.ap0_fft_v(h, z.ctypes.data, Z.ctypes.data, npair)
    v = Z.view(np.float64).reshape(npair, N, 2)
    P = np.einsum('ijk,ijk->j', v, v)
    return finish(P, fr, m)

def welch_ctrl(x):   # same glue, pocketfft complex FFT (scipy.fft, 1 worker)
    fr, m = glue_pack(x)
    Zc = scipy.fft.fft(z, axis=-1, workers=1)
    v = Zc.view(np.float64).reshape(npair, N, 2)
    P = np.einsum('ijk,ijk->j', v, v)
    return finish(P, fr, m)

def welch_ref(x):
    return scipy.signal.welch(x, fs=fs, window='hann', nperseg=N)[1]

def welch_ref4(x):
    with scipy.fft.set_workers(4):
        return scipy.signal.welch(x, fs=fs, window='hann', nperseg=N)[1]

rng = np.random.default_rng(7)
x = rng.standard_normal(L) + 0.3 * np.sin(0.05 * np.arange(L))
r = welch_ref(x); a = welch_sme(x); c = welch_ctrl(x)
ea = np.max(np.abs(a - r)) / np.max(np.abs(r)); ec = np.max(np.abs(c - r)) / np.max(np.abs(r))
print(f"numpy {np.__version__} scipy {scipy.__version__}; N={N} L=2^{sys.argv[2]} nseg={nseg}")
print(f"CHECK sme vs scipy.welch max rel {ea:.2e} {'PASS' if ea < 1e-12 else 'FAIL'}; ctrl {ec:.2e}")
# broken copy: wrong pairing sign must fail
saved = neg.copy(); neg[:] = np.arange(N); b = welch_sme(x); neg[:] = saved
eb = np.max(np.abs(b - r)) / np.max(np.abs(r))
print(f"CHECK broken (no conj-mirror fold) max rel {eb:.2e} {'FAIL as required' if eb > 1e-6 else 'NOT CAUGHT'}")
if mode != "time": sys.exit(0)

progs = [("scipy.signal.welch (1 thread)", lambda: welch_ref(x)), ("scipy.signal.welch set_workers(4)", lambda: welch_ref4(x)),
         ("same glue + pocketfft c2c (control)", lambda: welch_ctrl(x)), ("same glue + SME fp64 (ap-0)", lambda: welch_sme(x))]
calls = {}
for name, f in progs:
    f(); k = 1
    while True:
        t = time.perf_counter()
        for _ in range(k): f()
        dt = time.perf_counter() - t
        if dt >= 0.2: break
        k = max(k * 2, int(k * 0.25 / max(dt, 1e-6)))
    calls[name] = k
reps = {n: [] for n, _ in progs}
for rr in range(9):
    for j in range(len(progs)):
        name, f = progs[(j + rr) % len(progs)]
        k = calls[name]; t = time.perf_counter()
        for _ in range(k): f()
        reps[name].append((time.perf_counter() - t) / k * 1e3)
base = sorted(reps[progs[0][0]])[4]
for name, _ in progs:
    v = sorted(reps[name])
    print(f"RESULT N={N} L=2^{sys.argv[2]} {name:40s} best {v[0]:8.3f} median {v[4]:8.3f} ms/call  vs scipy median {100*(v[4]/base-1):+6.1f} %  ns/sample {v[4]*1e6/L:6.3f}")
