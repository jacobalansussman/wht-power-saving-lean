# as0py.py -- ctypes wrapper of libas0.dylib (wave 14, as-0): float64 Welch PSD with the SME, vDSP or PFFFT-d engine.
import ctypes, os
import numpy as np, scipy.signal
_lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "libas0.dylib"))
_lib.as0_new.restype = ctypes.c_void_p
_lib.as0_new.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_void_p, ctypes.c_double]
_lib.as0_welch_rows.restype = ctypes.c_int
_lib.as0_welch_rows.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_void_p, ctypes.c_int]
ENGINES = {"sme": 0, "smef": 1, "vdsp": 2, "pfd": 3, "smed": 4, "kfr": 5, "smew": 6, "fftw": 7}
_cache = {}

def welch(x, fs=1.0, window="hann", nperseg=256, noverlap=None, engine="sme"):
    """scipy.signal.welch(x, fs, window, nperseg, noverlap) for float64 x (1-D or 2-D rows, last axis),
    detrend='constant', scaling='density', average='mean', return_onesided=True."""
    if noverlap is None: noverlap = nperseg // 2
    key = (nperseg, noverlap, window if isinstance(window, str) else id(window), fs)
    h = _cache.get(key)
    if h is None:
        win = np.ascontiguousarray(scipy.signal.get_window(window, nperseg), dtype=np.float64)
        h = _lib.as0_new(nperseg, noverlap, win.ctypes.data, fs); _cache[key] = h
    x = np.ascontiguousarray(x, dtype=np.float64)
    shp = x.shape; X = x.reshape(-1, shp[-1])
    out = np.empty((X.shape[0], nperseg // 2 + 1))
    rc = _lib.as0_welch_rows(h, X.ctypes.data, X.shape[0], X.shape[1], out.ctypes.data, ENGINES[engine])
    if rc: raise ValueError(f"as0 engine {engine} rc={rc} (unsupported size)")
    f = np.fft.rfftfreq(nperseg, 1.0 / fs)
    return f, out.reshape(shp[:-1] + (nperseg // 2 + 1,))
