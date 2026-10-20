"""fc64.py -- av-1 (wave 15): drop-in np.convolve / np.correlate for float64 and complex128 that run the direct
sum on the M4's SME2 matrix unit (libfc64.dylib, Toeplitz tiles), NEON for short calls, and an FFT above a
crossover (scipy.signal.oaconvolve). Same modes and argument rules as numpy (full / same / valid; numpy swaps
the arguments of convolve when the kernel is longer than the signal; correlate conjugates its second argument).

Primitive (C): y[t] = sum_k g[k] * xp[o0 + t + k], xp = x padded with L-1 zeros on the left (and zeros right);
t indexes numpy's FULL output.
"""
import ctypes
import os

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.environ.get("FC64_LIB", os.path.join(_here, "libfc64.dylib")))
_P = ctypes.c_void_p
_L = ctypes.c_long
for _f in (_lib.fc_corr_d, _lib.fc_corr_z):
    _f.argtypes = [_P, _L, _P, ctypes.c_int, _L, _L, _P, ctypes.c_int]
    _f.restype = ctypes.c_int
_lib.fc_vdsp_d.argtypes = [_P, _L, _P, ctypes.c_int, _L, _L, _P]
_lib.fc_vdsp_d.restype = ctypes.c_int
_corr_d = _lib.fc_corr_d
_corr_z = _lib.fc_corr_z
_vdsp_d = _lib.fc_vdsp_d

PATH = 0  # 0 auto, 1 SME, 2 NEON, 3 scalar
FFT_MIN_TAPS = 2400  # float64: scipy oaconvolve at or above this many taps (t8: crossover against scipy's FFT routes ~2,350)
FFT_MIN_TAPS_C = 1100  # complex128: crossover ~1,100 taps (t7, t9: -7 / -8 % at 1,025)


def _glib():
    """the loaded image whose globals the active binding uses (the extension's own copy when it is in use)"""
    return ctypes.CDLL(_fc64ext.__file__) if _ext is not None else _lib


def set_global(name, value):
    ctypes.c_long.in_dll(_glib(), name).value = value


def set_broken(flag):
    ctypes.c_int.in_dll(_glib(), "fc_broken").value = 1 if flag else 0


def _range(n, m, mode):
    """(o0, nout) inside the full output of length n + m - 1, numpy's rule when n >= m (no swap)."""
    if mode == "full":
        return 0, n + m - 1
    if mode == "same":
        return (m - 1) // 2, n
    if mode == "valid":
        return m - 1, n - m + 1
    raise ValueError(mode)


try:  # the CPython extension (same C code, buffer protocol, GIL released); ctypes is the fallback
    if os.environ.get("FC64_CTYPES"):
        raise ImportError
    import sys as _sys
    _sys.path.insert(0, os.environ.get("FC64_EXTDIR", _here))
    import _fc64ext
    _ext = _fc64ext.corr
except ImportError:
    _ext = None
_F8 = np.dtype(np.float64)
_C16 = np.dtype(np.complex128)
_ND = np.ndarray


def _run(x, g, o0, nout, cplx, fn=None, path=None):
    """y[t] = sum_k g[k] xp[o0+t+k]; x, g contiguous float64 / complex128."""
    path = PATH if path is None else path
    y = np.empty(nout, dtype=_C16 if cplx else _F8)
    if nout <= 0:
        return y
    if _ext is not None:
        _ext(x, g, o0, y, 9 if fn is _vdsp_d else path)
        return y
    if fn is None:
        fn = _corr_z if cplx else _corr_d
    if fn is _vdsp_d:
        r = fn(x.ctypes.data, x.shape[0], g.ctypes.data, g.shape[0], o0, nout, y.ctypes.data)
    else:
        r = fn(x.ctypes.data, x.shape[0], g.ctypes.data, g.shape[0], o0, nout, y.ctypes.data, path)
    if r:
        raise MemoryError("fc64 scratch allocation failed")
    return y


def _prep(a, v):
    a = np.asarray(a)
    v = np.asarray(v)
    if a.ndim != 1 or v.ndim != 1:
        raise ValueError("object too deep for desired array")
    if a.size == 0:
        raise ValueError("a cannot be empty")
    if v.size == 0:
        raise ValueError("v cannot be empty")
    rt = np.result_type(a, v)
    if rt == np.float64:
        cplx = False
    elif rt == np.complex128:
        cplx = True
    else:
        return None  # other dtypes: numpy's own routine
    dt = rt
    return np.ascontiguousarray(a, dtype=dt), np.ascontiguousarray(v, dtype=dt), cplx


def convolve(a, v, mode="full", fn=None, path=None):
    # fast path: two 1-D contiguous float64 arrays (the common call)
    if (type(a) is _ND and type(v) is _ND and a.dtype is _F8 and v.dtype is _F8 and a.ndim == 1 and v.ndim == 1
            and a.flags.c_contiguous and v.flags.c_contiguous and a.shape[0] and v.shape[0]):
        cplx = False
    else:
        pr = _prep(a, v)
        if pr is None:
            return np.convolve(a, v, mode)
        a, v, cplx = pr
    if v.shape[0] > a.shape[0]:
        a, v = v, a
    n, m = a.shape[0], v.shape[0]
    o0, nout = _range(n, m, mode)
    if fn is None and path is None:  # the measured dispatch (av-1 wave 15, one P-core, whole call)
        d = _choose(nout, m, cplx)
        if d == "np":
            return np.convolve(a, v, mode)
        if d == "fft":
            import scipy.signal as ss
            return ss.oaconvolve(a, v, mode)
        if d == "vdsp":
            fn = _vdsp_d
        path = 1 if cplx else 5  # float64: the blocked kernel (error at numpy's level, +0-7 % time, t12/t13)
    g = np.ascontiguousarray(v[::-1])
    return _run(a, g, o0, nout, cplx, fn, path)


def _choose(nout, m, cplx):
    """Route per (outputs, taps), from the short-call and grid timings (t2, t5, t6, t8 files):
    SME direct where it measured fastest; numpy's own loop for tiny kernels on short calls; vDSP_convD for short calls
    with longer kernels (it beat both SME and our NEON loop at 256 outputs); scipy's overlap-add FFT above FFT_MIN_TAPS."""
    if m >= (FFT_MIN_TAPS_C if cplx else FFT_MIN_TAPS):
        return "fft"
    if cplx:
        return "sme" if (nout >= 4096 and m >= 9) or nout * m >= 300000 else "np"
    if (m >= 9 and nout >= 4096) or (m >= 33 and nout >= 1024) or (m < 9 and nout >= 100000):
        return "sme"
    if m < 9 or (m < 33 and nout < 1024 and m <= 9):
        return "np"
    return "vdsp"


def correlate(a, v, mode="valid", fn=None, path=None):
    pr = _prep(a, v)
    if pr is None:
        return np.correlate(a, v, mode)
    a, v, cplx = pr
    n, m = a.shape[0], v.shape[0]
    if m <= n:
        o0, nout = _range(n, m, mode)
        g = v.conj() if cplx else v
        return _run(a, g, o0, nout, cplx, fn, path)
    # kernel longer than the signal: numpy's rule (inverted) for the lag range, computed as conj(reverse(c(v, a)))
    if mode == "full":
        o0, nout = 0, n + m - 1
    elif mode == "same":
        o0, nout = n // 2, m
    else:
        o0, nout = n - 1, m - n + 1
    s_lo = (n + m - 2) - (o0 + nout - 1)
    g = a.conj() if cplx else a
    ys = _run(v, g, s_lo, nout, cplx, fn, path)
    ys = ys[::-1]
    return ys.conj() if cplx else ys.copy()


def vdsp_convolve(a, v, mode="full"):
    """Baseline: Accelerate vDSP_convD through the same wrapper (float64 only)."""
    return convolve(a, v, mode, fn=_vdsp_d)
