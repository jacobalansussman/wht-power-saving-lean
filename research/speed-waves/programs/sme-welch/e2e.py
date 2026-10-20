# e2e.py -- whole calls from Python, float64, C channels x L samples (C-contiguous 2-D array):
#   scipy.signal.welch (1 thread; it has no workers argument, scipy.fft.set_workers(4) is the only threading knob),
#   mne.time_frequency.psd_array_welch as users call it (n_jobs=1 and n_jobs=4),
#   libas0 through ctypes: vDSP-D pipeline (C baseline) and the SME engine (SME-F at 512/4096, SME-D at 1024/2048).
# usage: python e2e.py N NOVERLAP C LOG2L WINDOW
import sys, os, time, numpy as np, scipy, scipy.signal, scipy.fft
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import as0py, mne
mne.set_log_level("ERROR")
N, nov, C, LL, win = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
L = 1 << LL; fs = 1000.0
sme = "smef" if N in (512, 4096) else "smed"
rng = np.random.default_rng(3)
x = rng.standard_normal((C, L)) + 0.2
def sc1(): return scipy.signal.welch(x, fs=fs, window=win, nperseg=N, noverlap=nov)[1]
def sc4():
    with scipy.fft.set_workers(4): return scipy.signal.welch(x, fs=fs, window=win, nperseg=N, noverlap=nov)[1]
def mne1(): return mne.time_frequency.psd_array_welch(x, fs, n_fft=N, n_overlap=nov, n_per_seg=N, window=win, average="mean", remove_dc=True, n_jobs=1)[0]
def mne4(): return mne.time_frequency.psd_array_welch(x, fs, n_fft=N, n_overlap=nov, n_per_seg=N, window=win, average="mean", remove_dc=True, n_jobs=4)[0]
def vd(): return as0py.welch(x, fs=fs, window=win, nperseg=N, noverlap=nov, engine="vdsp")[1]
def sm(): return as0py.welch(x, fs=fs, window=win, nperseg=N, noverlap=nov, engine=sme)[1]
r = sc1()
print(f"numpy {np.__version__} scipy {scipy.__version__} mne {mne.__version__}; N={N} nov={nov} C={C} L=2^{LL} window={win}")
for name, f in (("scipy w4", sc4), ("mne n_jobs=1", mne1), ("as0 vdsp", vd), ("as0 " + sme, sm)):
    a = f(); e = np.max(np.abs(a - r) / np.abs(r))
    print(f"CHECK {name:14s} vs scipy.signal.welch: max bin rel {e:.2e} {'PASS' if e < 1e-13 else 'FAIL'}")
if os.environ.get("CHECK_ONLY"): sys.exit(0)
progs = [("scipy.signal.welch 1 thread", sc1), ("scipy.signal.welch set_workers(4)", sc4), ("MNE psd_array_welch n_jobs=1", mne1),
         ("MNE psd_array_welch n_jobs=4", mne4), ("C vDSP-D pipeline (ctypes)", vd), (f"SME {sme} (ctypes)", sm)]
if os.environ.get("NO_MNE4"): progs = [p for p in progs if "n_jobs=4" not in p[0]]
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
med = {n: sorted(v)[4] for n, v in reps.items()}
base = med[progs[0][0]]; vdm = med["C vDSP-D pipeline (ctypes)"]
for name, _ in progs:
    v = sorted(reps[name])
    print(f"RESULT N={N} nov={nov} C={C} L=2^{LL} {win} {name:36s} best {v[0]:9.3f} median {v[4]:9.3f} ms/call  ns/sample {v[4]*1e6/(C*L):6.3f}"
          f"  vs scipy-1 {100*(v[4]/base-1):+6.1f} %  vs C-vDSP {100*(v[4]/vdm-1):+6.1f} %")
