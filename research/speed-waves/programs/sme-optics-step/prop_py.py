# prop_py.py -- as-1: angular-spectrum propagation step, end to end from Python (numpy 2.0.2, scipy 1.13.1).
# Workload as scout-c0's pipes2.py: E random complex128 N x N, H = exp(i 2 pi z sqrt(1/lam^2 - fx^2 - fy^2)),
# dx = 10 um, lam = 633 nm, z = 0.1 m; asp() = ifft2(fft2(E) * H).
# usage: python prop_py.py check N | python prop_py.py time N [progs...]
# Timing: >= 0.2 s per repeat, 9 repeats, programs take turns (rotating start), ms per propagation step.
import os, sys, time, ctypes
import numpy as np, scipy.fft as sf

N = int(sys.argv[2]); mode = sys.argv[1]
rng = np.random.default_rng(2)
E = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
fx = np.fft.fftfreq(N, 1e-5); FX, FY = np.meshgrid(fx, fx, indexing="ij"); lam = 633e-9
H = np.exp(1j * 2 * np.pi * 0.1 * np.sqrt(np.maximum(0, 1 / lam ** 2 - FX ** 2 - FY ** 2)))
del FX, FY
lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "libas1.dylib"))
lib.as1_new.restype = ctypes.c_void_p; lib.as1_new.argtypes = [ctypes.c_long, ctypes.c_void_p]
for f in (lib.as1_prop, lib.as1_prop_nat, lib.as1_fft2): f.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
lib.as1_step.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
if os.environ.get("PFD"): lib.as1_set_pfd(int(os.environ["PFD"]))
if os.environ.get("NOFUSE"): lib.as1_set_fused(0)
h = lib.as1_new(N, H.ctypes.data)
P = lambda a: a.ctypes.data

def sme_prop_rev():  # fresh output, reversed order (exact in chains of even length)
    out = np.empty_like(E); lib.as1_prop(h, P(E), P(out)); return out
def sme_prop_nat():  # fresh output, natural order (one extra reversal pass)
    out = np.empty_like(E); lib.as1_prop_nat(h, P(E), P(out)); return out
U = E.copy()
def sme_step_inplace(): lib.as1_step(h, P(U)); return U
def np_asp(): return np.fft.ifft2(np.fft.fft2(E) * H)
def sp1(): return sf.ifft2(sf.fft2(E, workers=1) * H, workers=1)
def sp4(): return sf.ifft2(sf.fft2(E, workers=4) * H, workers=4)
def sp4o():  # scipy, 4 workers, intermediate reused (multiply in place, overwrite_x)
    a = sf.fft2(E, workers=4); a *= H; return sf.ifft2(a, workers=4, overwrite_x=True)
def sp1o():
    a = sf.fft2(E, workers=1); a *= H; return sf.ifft2(a, workers=1, overwrite_x=True)
lib.as1_mt_init.argtypes = [ctypes.c_void_p]; lib.as1_prop_mt.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int]
if os.environ.get("MT"): lib.as1_mt_init(h)
def sme_mt():  # 1 SME + 3 NEON threads, fresh output, reversed order
    out = np.empty_like(E); lib.as1_prop_mt(h, P(E), P(out), 1); return out
def neon_mt():  # 4 NEON threads, same plan
    out = np.empty_like(E); lib.as1_prop_mt(h, P(E), P(out), 0); return out
progs = {"cluster 1 SME + 3 NEON (rev.)": sme_mt, "cluster 4 NEON (rev.)": neon_mt, "SME prop (rev.)": sme_prop_rev, "SME prop (natural)": sme_prop_nat, "SME step in place (rev.)": sme_step_inplace,
         "numpy asp": np_asp, "scipy w=1": sp1, "scipy w=1 overwrite": sp1o, "scipy w=4": sp4, "scipy w=4 overwrite": sp4o}

if mode == "check":
    ref = np_asp(); nrm = np.abs(ref).max()
    r = sme_prop_nat(); print(f"CHECK SME prop natural vs numpy asp: max rel {np.abs(r - ref).max() / nrm:.3e}  rel RMS {np.linalg.norm(r - ref) / np.linalg.norm(ref):.3e}")
    r = sme_prop_rev(); rr = np.roll(r[::-1, ::-1], 1, axis=(0, 1))
    print(f"CHECK SME prop reversed (read at -n) vs numpy asp: max rel {np.abs(rr - ref).max() / nrm:.3e}")
    if os.environ.get("MT"):
        for f, nm in ((sme_mt, "1 SME + 3 NEON"), (neon_mt, "4 NEON")):
            r = f(); rr = np.roll(r[::-1, ::-1], 1, axis=(0, 1)); print(f"CHECK cluster {nm} (read at -n) vs numpy asp: max rel {np.abs(rr - ref).max() / nrm:.3e}")
    r2 = sp4(); print(f"CHECK scipy w=4 vs numpy asp: max rel {np.abs(r2 - ref).max() / nrm:.3e}")
    y = np.empty_like(E); lib.as1_fft2(h, P(E), P(y)); f = np.fft.fft2(E)
    print(f"CHECK SME fft2 vs numpy fft2: max rel {np.abs(y - f).max() / np.abs(f).max():.3e}  rel RMS {np.linalg.norm(y - f) / np.linalg.norm(f):.3e}")
    sys.exit(0)

names = sys.argv[3:] or [n for n in progs if os.environ.get("MT") or not n.startswith("cluster")]
print(f"# numpy {np.__version__} scipy {sf.__name__} {__import__('scipy').__version__}; N={N}; ms per propagation step")
calls = {}
for n in names:
    f = progs[n]; f(); k = 1
    while True:
        t0 = time.perf_counter(); [f() for _ in range(k)]; dt = time.perf_counter() - t0
        if dt >= 0.2: break
        k = max(k + 1, int(k * 0.22 / max(dt, 1e-6)) + 1)
    calls[n] = k
res = {n: [] for n in names}
for r in range(9):
    for jj in range(len(names)):
        n = names[(jj + r) % len(names)]; f = progs[n]; k = calls[n]
        t0 = time.perf_counter(); [f() for _ in range(k)]; dt = (time.perf_counter() - t0) / k * 1e3
        res[n].append(dt); print(f"rep {r} {n:28s} {dt:9.3f}")
for n in names:
    v = sorted(res[n]); print(f"RESULT N={N} {n:28s} best {v[0]:9.3f}  median {v[4]:9.3f}  ms/step  (calls/rep {calls[n]})")
