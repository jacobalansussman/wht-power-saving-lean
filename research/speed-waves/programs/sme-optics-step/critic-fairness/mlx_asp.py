# au-f: angular-spectrum step ifft2(fft2(E)*H) on the M4 GPU with MLX, complex64 (the GPU has no fp64).
# Relevance check for as-1 (complex128 on SME): what a user who accepts fp32 gets on this Mac.
import time, sys, numpy as np, mlx.core as mx
def H_of(N, dx=10e-6, lam=633e-9, z=0.1):
    f = np.fft.fftfreq(N, dx); fx, fy = np.meshgrid(f, f, indexing='xy')
    arg = 1.0/lam**2 - fx**2 - fy**2
    kz = 2*np.pi*np.sqrt(np.maximum(arg, 0.0)); H = np.exp(1j*kz*z) * (arg > 0)
    return H
rng = np.random.default_rng(1)
for N in [int(a) for a in sys.argv[1:]] or [1024, 2048, 4096]:
    E = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))
    H = H_of(N)
    Eg = mx.array(E.astype(np.complex64)); Hg = mx.array(H.astype(np.complex64)); mx.eval(Eg, Hg)
    step = lambda: mx.fft.ifft2(mx.fft.fft2(Eg) * Hg)
    out = step(); mx.eval(out)
    if N <= 2048:
        ref = np.fft.ifft2(np.fft.fft2(E) * H)
        o = np.array(out).astype(np.complex128)
        print(f"CHECK N={N} GPU c64 vs numpy c128: max rel {np.abs(o-ref).max()/np.abs(ref).max():.3e}  rel RMS {np.sqrt((np.abs(o-ref)**2).sum()/(np.abs(ref)**2).sum()):.3e}")
        # broken copy: conj of H must be far off
        ob = np.array(mx.fft.ifft2(mx.fft.fft2(Eg) * mx.conj(Hg))).astype(np.complex128)
        print(f"CHECK N={N} broken (conj H): rel RMS {np.sqrt((np.abs(ob-ref)**2).sum()/(np.abs(ref)**2).sum()):.3e} (must be large)")
    # timing: 9 repeats of >= 0.2 s each
    k = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(k): r = step(); mx.eval(r)
        if time.perf_counter() - t0 >= 0.2: break
        k *= 2
    ts = []
    for rep in range(9):
        t0 = time.perf_counter()
        for _ in range(k): r = step(); mx.eval(r)
        ts.append((time.perf_counter() - t0) / k * 1e3)
        print(f"rep {rep} N={N} GPU c64 step {ts[-1]:.3f} ms")
    ts.sort(); print(f"RESULT N={N} MLX GPU complex64 step best {ts[0]:.3f} median {ts[4]:.3f} ms (calls/rep {k})", flush=True)
    del Eg, Hg, out
