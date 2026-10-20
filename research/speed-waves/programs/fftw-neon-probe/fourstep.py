# au-h0 probe: numpy's own FFT, one long 1-D complex128 transform, against a four-step built from numpy's own
# batched calls (pocketfft vectorises across rows/columns in batched calls; a single 1-D call is not vectorised).
import time, numpy as np, sys
print("numpy", np.__version__)
rng = np.random.default_rng(1)
def best(f, reps=9):
    ts = []
    for _ in range(reps):
        it = 1
        while True:
            t0 = time.perf_counter()
            for _ in range(it): f()
            t = time.perf_counter() - t0
            if t > 0.2: break
            it *= 2
        ts.append(t / it)
    ts.sort(); return ts[0], ts[len(ts)//2]
for lg in [int(a) for a in sys.argv[1:]]:
    n = 1 << lg; n1 = 1 << (lg // 2); n2 = n // n1
    x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    k1 = np.arange(n1)[:, None]; j2 = np.arange(n2)[None, :]
    tw = np.exp(-2j * np.pi * k1 * j2 / n)            # n1 x n2
    def four():
        a = x.reshape(n2, n1).T                        # a[j1, j2] = x[j1 + n1*j2]... use x[n1*j2 + j1]
        b = np.fft.fft(a, axis=1)                      # length-n2 transforms over j2 -> k2 (batched, strided)
        b *= tw                                        # twiddle w^(j1*k2)
        c = np.fft.fft(b, axis=0)                      # length-n1 transforms over j1 -> k1
        return c.reshape(-1)                           # X[k2 + n2*k1] = c[k1, k2]
    ref = np.fft.fft(x)
    y = four()
    err = np.max(np.abs(y - ref)) / np.max(np.abs(ref))
    tw = np.exp(-2j * np.pi * np.arange(n1)[:, None] * np.arange(n2)[None, :] / n)
    a1 = best(lambda: np.fft.fft(x)); a2 = best(four)
    print(f"n=2^{lg} ({n1}x{n2}) err {err:.2e} | np.fft.fft single: best {a1[0]/n*1e9:.3f} med {a1[1]/n*1e9:.3f} ns/pt"
          f" | four-step from numpy's batched calls: best {a2[0]/n*1e9:.3f} med {a2[1]/n*1e9:.3f} ns/pt | ratio {a2[0]/a1[0]:.2f}")
