# ap0_backend.py -- ap-0: SKETCH of a scipy.fft backend (uarray protocol, as pyFFTW and mkl_fft do) that sends
# forward complex128 FFTs along the last axis of a C-contiguous array, n = 512, 1024, 2048, 4096 (no padding/cropping),
# norm=None/"backward", to the SME fp64 kernel; everything else returns NotImplemented, so scipy falls back to pocketfft.
# Use:  import scipy.fft, ap0_backend
#       with scipy.fft.set_backend(ap0_backend, coerce=False, only=False): y = scipy.fft.fft(x)
import ctypes, os
import numpy as np
_lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "libap0.dylib"))
_lib.ap0_new_v.restype = ctypes.c_void_p; _lib.ap0_new_v.argtypes = [ctypes.c_long]
_lib.ap0_fft_v.restype = None; _lib.ap0_fft_v.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long]
_plans = {}
SIZES = (512, 1024, 2048, 4096)
__ua_domain__ = "numpy.scipy.fft"
calls = {"sme": 0, "fallback": 0}

def _fft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *, plan=None):
    x = np.asarray(x)
    if (x.dtype != np.complex128 or x.ndim < 1 or (n is not None and n != x.shape[axis]) or axis not in (-1, x.ndim - 1)
            or norm not in (None, "backward") or not x.flags.c_contiguous or x.shape[-1] not in SIZES or plan is not None):
        calls["fallback"] += 1
        return NotImplemented
    N = x.shape[-1]; rows = x.size // N
    h = _plans.get(N) or _plans.setdefault(N, _lib.ap0_new_v(N))
    out = np.empty_like(x)  # numpy semantics: a new array (allocation counted)
    _lib.ap0_fft_v(h, x.ctypes.data, out.ctypes.data, rows)
    calls["sme"] += 1
    return out

def __ua_function__(method, args, kwargs):
    if method.__name__ == "fft":
        return _fft(*args, **kwargs)
    calls["fallback"] += 1
    return NotImplemented
