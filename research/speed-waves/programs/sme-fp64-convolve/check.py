"""check.py -- correctness of fc64 against numpy and an exact (correctly rounded) reference on sampled outputs.
Usage: python check.py [quick|full|broken]
Exact reference: every product g[k]*x[i] split exactly into p + e (Dekker / Veltkamp, no overflow at these
magnitudes), then math.fsum over all 2L (complex: 4L per part) terms = the correctly rounded sum.
Error per output is reported as |got - exact| / sum_k |g[k] x[i]| (the condition-free scale of a direct sum),
max over the sampled outputs; numpy's own error is computed on the same outputs.
"""
import math
import sys

import numpy as np

import fc64

rng = np.random.default_rng(12345)
SPLIT = 134217729.0  # 2^27 + 1


def two_prod(a, b):
    p = a * b
    ca = SPLIT * a
    ah = ca - (ca - a)
    al = a - ah
    cb = SPLIT * b
    bh = cb - (cb - b)
    bl = b - bh
    e = ((ah * bh - p) + ah * bl + al * bh) + al * bl
    return p, e


def exact_full(x, g, ts):
    """full-convolution-index outputs y[t] = sum_k g[k] xp[t+k] (correlation form, xp = x padded by L-1)."""
    L = len(g)
    n = len(x)
    out, scale = [], []
    for t in ts:
        ks = np.arange(L)
        ii = t + ks - (L - 1)
        ok = (ii >= 0) & (ii < n)
        gk = g[ks[ok]]
        xi = x[ii[ok]]
        if np.iscomplexobj(x):
            terms_r, terms_i = [], []
            for (a, b, sgn, lst) in ((gk.real, xi.real, 1, terms_r), (gk.imag, xi.imag, -1, terms_r),
                                     (gk.real, xi.imag, 1, terms_i), (gk.imag, xi.real, 1, terms_i)):
                p, e = two_prod(a, b)
                lst.extend((sgn * p).tolist())
                lst.extend((sgn * e).tolist())
            out.append(complex(math.fsum(terms_r), math.fsum(terms_i)))
            scale.append(float(np.sum(np.abs(gk) * np.abs(xi))))
        else:
            p, e = two_prod(gk, xi)
            out.append(math.fsum(p.tolist() + e.tolist()))
            scale.append(float(np.sum(np.abs(gk * xi))))
    return np.array(out), np.array(scale)


def sample_ts(nout, k=120):
    ts = list(range(min(nout, 16))) + list(range(max(0, nout - 16), nout))
    if nout > 32:
        ts += rng.integers(0, nout, size=k).tolist()
    return sorted(set(ts))


def one(kind, n, m, mode, cplx, path, results):
    if cplx:
        a = rng.uniform(-1, 1, n) + 1j * rng.uniform(-1, 1, n)
        v = rng.uniform(-1, 1, m) + 1j * rng.uniform(-1, 1, m)
    else:
        a = rng.uniform(-1, 1, n)
        v = rng.uniform(-1, 1, m)
    if kind == "convolve":
        ref_np = np.convolve(a, v, mode)
        got = fc64.convolve(a, v, mode, path=path)
    else:
        ref_np = np.correlate(a, v, mode)
        got = fc64.correlate(a, v, mode, path=path)
    if got.shape != ref_np.shape or got.dtype != ref_np.dtype:
        results.append(False)
        print(f"FAIL shape/dtype {kind} n={n} m={m} {mode} cplx={cplx} path={path}: {got.shape} {got.dtype} vs {ref_np.shape} {ref_np.dtype}")
        return
    d = np.max(np.abs(got - ref_np)) / max(np.max(np.abs(ref_np)), 1e-300)
    # exact reference on sampled outputs: map numpy's output index to the full correlation-form index
    if kind == "convolve":
        x, w = (a, v) if n >= m else (v, a)
        g = w[::-1]
        L = len(w)
        o0 = fc64._range(len(x), L, mode)[0]
        ts = sample_ts(len(got))
        ex, sc = exact_full(x, g, [o0 + t for t in ts])
    else:  # correlate: c[k] = sum_n a[n+k] conj(v[n]); full index t <-> lag t-(m-1)
        g = v.conj()
        if n >= m:
            o0 = fc64._range(n, m, mode)[0]
        else:
            o0 = {"full": 0, "same": n // 2, "valid": n - 1}[mode]
        ts = sample_ts(len(got))
        ex, sc = exact_full(a, g, [o0 + t for t in ts])
    sc = np.maximum(sc, 1e-300)
    e_got = np.max(np.abs(got[ts] - ex) / sc)
    e_np = np.max(np.abs(ref_np[ts] - ex) / sc)
    ok = d < 1e-12 and e_got < 1e-13
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'} {kind:9s} n={n:7d} m={m:5d} {mode:5s} {'c128' if cplx else 'f64 '} path={path} "
          f"vs numpy {d:.2e}  err/exact: fc64 {e_got:.2e} numpy {e_np:.2e}  ratio {e_got / max(e_np, 1e-300):.2f}")
    return e_got, e_np


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "quick"
    res = []
    if which == "broken":
        fc64.set_broken(True)
        for (n, m) in ((5000, 33), (100000, 257)):
            one("convolve", n, m, "same", False, 1, res)
            one("convolve", n, m, "same", True, 1, res)
        fc64.set_broken(False)
        print(f"BROKEN COPY: {sum(not r for r in res)} of {len(res)} checks FAIL (all must FAIL)")
        return
    sizes = [(1, 1), (3, 1), (7, 3), (100, 3), (100, 4), (513, 33), (1000, 1025), (4096, 64), (5000, 33),
             (10007, 128), (20000, 257), (65536, 1025), (12345, 1024), (999, 1000)]
    if which == "full":
        sizes += [(100000, 3), (100000, 9), (100000, 513), (300001, 129), (1000000, 33), (200000, 2049)]
    for (n, m) in sizes:
        for mode in ("full", "same", "valid"):
            for cplx in (False, True):
                for kind in ("convolve", "correlate"):
                    paths = (1, 2, 4, 5) if not cplx else (1,)
                    for path in paths:
                        if path == 2 and n * m > 3e7:
                            continue
                        one(kind, n, m, mode, cplx, path, res)
                        if n != m:  # swapped arguments (kernel longer than the signal)
                            pass
        # swapped: kernel longer than the signal
        for mode in ("full", "same", "valid"):
            for cplx in (False, True):
                for kind in ("convolve", "correlate"):
                    if m < n and n * m < 3e7:
                        rng_state = None
                        one_swapped(kind, n, m, mode, cplx, res)
    print(f"SUMMARY {sum(res)} of {len(res)} PASS")


def one_swapped(kind, n, m, mode, cplx, res):
    # a has m samples, v has n (> m): numpy swaps
    if cplx:
        a = rng.uniform(-1, 1, m) + 1j * rng.uniform(-1, 1, m)
        v = rng.uniform(-1, 1, n) + 1j * rng.uniform(-1, 1, n)
    else:
        a = rng.uniform(-1, 1, m)
        v = rng.uniform(-1, 1, n)
    f_np = np.convolve if kind == "convolve" else np.correlate
    f_fc = fc64.convolve if kind == "convolve" else fc64.correlate
    r = f_np(a, v, mode)
    g = f_fc(a, v, mode, path=1)
    ok = g.shape == r.shape and g.dtype == r.dtype and np.max(np.abs(g - r)) / max(np.max(np.abs(r)), 1e-300) < 1e-12
    res.append(ok)
    d = np.max(np.abs(g - r)) / max(np.max(np.abs(r)), 1e-300) if g.shape == r.shape else float("nan")
    print(f"{'PASS' if ok else 'FAIL'} {kind:9s} SWAPPED len(a)={m} len(v)={n} {mode:5s} {'c128' if cplx else 'f64 '} vs numpy {d:.2e}")


if __name__ == "__main__":
    main()
