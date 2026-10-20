# e2e.py -- ap-0: end to end from Python (numpy 2.0.2, scipy 1.13.1): batched row FFTs of a 2-D complex128 array,
# the SME fp64 kernel through ctypes vs numpy.fft.fft and scipy.fft.fft (pocketfft). One thread unless stated.
# usage: python e2e.py N ROWS [N ROWS ...]
import ctypes, sys, time, os
import numpy as np, scipy, scipy.fft
lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "libap0.dylib"))
lib.ap0_new_v.restype = ctypes.c_void_p; lib.ap0_new_v.argtypes = [ctypes.c_long]
lib.ap0_fft_v.restype = None; lib.ap0_fft_v.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long]
print("numpy", np.__version__, "scipy", scipy.__version__, flush=True)

def bench(title, progs, npts):
    # calibrate: one repeat >= 0.2 s; 9 repeats, programs take turns
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
    reps = {name: [] for name, _ in progs}
    for r in range(9):
        for j in range(len(progs)):
            name, f = progs[(j + r) % len(progs)]
            k = calls[name]; t = time.perf_counter()
            for _ in range(k): f()
            reps[name].append((time.perf_counter() - t) / k / npts * 1e9)
    for name, _ in progs:
        v = sorted(reps[name])
        print(f"RESULT {title:40s} {name:44s} best {v[0]:8.4f} median {v[4]:8.4f} ns/point", flush=True)

args = sys.argv[1:]
for i in range(0, len(args), 2):
    N, R = int(args[i]), int(args[i + 1])
    rng = np.random.default_rng(1)
    x = (rng.standard_normal((R, N)) + 1j * rng.standard_normal((R, N))).astype(np.complex128)
    out = np.empty_like(x); y = np.empty_like(x)
    h = lib.ap0_new_v(N)
    xp, op = x.ctypes.data, out.ctypes.data
    def ours(): lib.ap0_fft_v(h, xp, op, R)
    def ours_alloc():
        o = np.empty_like(x); lib.ap0_fft_v(h, xp, o.ctypes.data, R); return o
    def np_fft(): return np.fft.fft(x, axis=-1)
    def sp1(): return scipy.fft.fft(x, axis=-1, workers=1)
    def sp4(): return scipy.fft.fft(x, axis=-1, workers=4)
    def sp1_ow():
        np.copyto(y, x); return scipy.fft.fft(y, axis=-1, workers=1, overwrite_x=True)
    ours(); ref = np.fft.fft(x, axis=-1)
    err = np.max(np.abs(out - ref)) / np.max(np.abs(ref))
    ow = sp1_ow(); same = ow.ctypes.data == y.ctypes.data
    print(f"CHECK N={N} rows={R}: SME vs numpy.fft max|diff|/max|X| = {err:.3e} {'PASS' if err < 1e-14 else 'FAIL'}; "
          f"scipy overwrite_x result in place: {same}", flush=True)
    progs = [("SME fp64 via ctypes, output reused", ours), ("SME fp64 via ctypes + np.empty_like per call", ours_alloc),
             ("numpy.fft.fft (allocates)", np_fft), ("scipy.fft.fft workers=1 (allocates)", sp1),
             ("scipy.fft.fft workers=1 copy+overwrite_x", sp1_ow), ("scipy.fft.fft workers=4 (allocates)", sp4)]
    bench(f"rows c128 N={N} rows={R} ({16*N*R/2**20:.0f} MiB)", progs, N * R)
