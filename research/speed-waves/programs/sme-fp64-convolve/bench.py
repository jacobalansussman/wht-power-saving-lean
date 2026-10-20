"""bench.py -- whole-Python-call timing of 1-D convolution routes on the M4 (av-1, wave 15).
Usage: python bench.py MODE DTYPE "N1,N2,..." "M1,M2,..." [progs]
  MODE = full|same|valid, DTYPE = f64|c128, progs = comma list (default: all that apply)
Programs (every call allocates its own output, like numpy):
  fc_sme    fc64.convolve, SME path (Toeplitz tiles)         fc_neon  fc64.convolve, NEON direct path
  fc_auto   fc64.convolve with the measured dispatch           np       numpy.convolve
  sp_auto   scipy.signal.convolve(method='auto')               sp_oa    scipy.signal.oaconvolve
  sp_fft    scipy.signal.fftconvolve                           vdsp     Accelerate vDSP_convD via the same ctypes wrapper
Timing (BRIEF section 5): perf_counter_ns (CLOCK_UPTIME_RAW on macOS), per program a call count so that one repeat
lasts >= 0.2 s, 9 repeats, programs take turns (the first program rotates each round), every repeat printed,
best and median in ns per output. Same input arrays again and again.
"""
import os
import sys
import time

import numpy as np
import scipy
import scipy.signal as ss

import fc64

MIN_S = 0.2
REPS = 9


def calib(f):
    f()
    k = 1
    while True:
        t0 = time.perf_counter_ns()
        for _ in range(k):
            f()
        dt = time.perf_counter_ns() - t0
        if dt >= MIN_S * 1e9:
            return k
        k = k * 8 if dt < MIN_S * 1e9 / 16 else int(np.ceil(k * MIN_S * 1e9 * 1.1 / dt))


def bench_cell(title, progs, nout):
    names = list(progs)
    ks = {nm: calib(progs[nm]) for nm in names}
    reps = {nm: [] for nm in names}
    for r in range(REPS):
        for jj in range(len(names)):
            nm = names[(jj + r) % len(names)]
            f, k = progs[nm], ks[nm]
            t0 = time.perf_counter_ns()
            for _ in range(k):
                f()
            dt = time.perf_counter_ns() - t0
            v = dt / (k * nout)
            reps[nm].append(v)
            print(f"rep {r:2d} {title} {nm:8s} {v:10.4f}")
    out = {}
    for nm in names:
        s = sorted(reps[nm])
        out[nm] = (s[0], s[len(s) // 2])
        print(f"RESULT {title} {nm:8s} best {s[0]:10.4f} median {s[len(s) // 2]:10.4f} ns/output (calls/rep {ks[nm]})")
    sys.stdout.flush()
    return out


_NC = None


def _noscan(f):
    """run f with the non-finite output scan off (fc_nancheck = 0), then switch it back on"""
    global _NC
    if _NC is None:
        import ctypes
        _NC = ctypes.c_int.in_dll(fc64._glib(), "fc_nancheck")
    _NC.value = 0
    try:
        return f()
    finally:
        _NC.value = 1


_FR = None


def _frac(fr, f):
    """run f with fc_frac (the SME share of the outputs in the multi-thread path) set to fr"""
    global _FR
    if _FR is None:
        import ctypes
        _FR = ctypes.c_double.in_dll(fc64._glib(), "fc_frac")
    _FR.value = fr
    return f()


def main():
    mode, dt = sys.argv[1], sys.argv[2]
    Ns = [int(float(s)) for s in sys.argv[3].split(",")]
    Ms = [int(float(s)) for s in sys.argv[4].split(",")]
    want = sys.argv[5].split(",") if len(sys.argv) > 5 else None
    print(f"# numpy {np.__version__} scipy {scipy.__version__} fc64 lib {fc64._lib._name} load {os.getloadavg()}")
    rng = np.random.default_rng(7)
    for N in Ns:
        for M in Ms:
            if dt == "c128":
                a = rng.uniform(-1, 1, N) + 1j * rng.uniform(-1, 1, N)
                v = rng.uniform(-1, 1, M) + 1j * rng.uniform(-1, 1, M)
            else:
                a = rng.uniform(-1, 1, N)
                v = rng.uniform(-1, 1, M)
            nout = len(np.convolve(a[: max(M, 8)], v, "full")) if False else {"full": N + M - 1, "same": max(N, M), "valid": max(N, M) - min(N, M) + 1}[mode]
            P = {
                "fc_sme": lambda: fc64.convolve(a, v, mode, path=1),
                "fc_neon": lambda: fc64.convolve(a, v, mode, path=2),
                "fc_smeb": lambda: fc64.convolve(a, v, mode, path=4),
                "fc_smeb8": lambda: fc64.convolve(a, v, mode, path=5),
                "fc_noscan": lambda: _noscan(lambda: fc64.convolve(a, v, mode, path=1)),
                "fc_mt3_75": lambda: _frac(0.75, lambda: fc64.convolve(a, v, mode, path=13)),
                "fc_mt3_85": lambda: _frac(0.85, lambda: fc64.convolve(a, v, mode, path=13)),
                "fc_mt1_85": lambda: _frac(0.85, lambda: fc64.convolve(a, v, mode, path=11)),
                "fc_b8noscan": lambda: _noscan(lambda: fc64.convolve(a, v, mode, path=5)),
                "fc_auto": lambda: fc64.convolve(a, v, mode),
                "np": lambda: np.convolve(a, v, mode),
                "sp_auto": lambda: ss.convolve(a, v, mode, method="auto"),
                "sp_oa": lambda: ss.oaconvolve(a, v, mode),
                "sp_fft": lambda: ss.fftconvolve(a, v, mode),
                "vdsp": lambda: fc64.vdsp_convolve(a, v, mode),
            }
            names = want or ["fc_sme", "np", "sp_auto", "sp_oa", "sp_fft", "vdsp"]
            if dt == "c128":
                names = [x for x in names if x not in ("vdsp", "fc_neon", "fc_smeb", "fc_smeb8")]
            progs = {nm: P[nm] for nm in names}
            # agreement check before timing
            ref = np.convolve(a, v, mode) if N * M <= 2e8 else ss.oaconvolve(a, v, mode)
            for nm in names:
                got = progs[nm]()
                err = np.max(np.abs(got - ref)) / np.max(np.abs(ref))
                print(f"CHECK N={N} M={M} {mode} {dt} {nm:8s} rel err vs {'np' if N * M <= 2e8 else 'oa'} {err:.2e} {'PASS' if err < 1e-12 else 'FAIL'}")
            bench_cell(f"N={N} M={M} {mode} {dt}", progs, nout)


if __name__ == "__main__":
    main()
