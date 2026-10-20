#!/usr/bin/env python3
"""dr3.py -- av-5: the same test as dr.py, also read with au-h1's (wave 14) stricter statistics: the MEAN over bins at least
40 bins from the tone (not the median), and the worst far bin against the float reference's worst far bin (spurs from
the Q14 twiddle table show here, not in the median)."""
import sys, numpy as np
sys.path.insert(0, '.')
from gen24 import model
rng = np.random.default_rng(11)
def run(L, K, fl):
    N = 1 << L; F = N // 2 + 1; t = np.arange(N * K); w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(N) / N)
    f0 = N / 12.3
    sig = 32767 * 10 ** (-6 / 20) * np.cos(2 * np.pi * f0 * t / N) + rng.normal(0, 32767 * 10 ** (fl / 20), N * K)
    q = np.clip(np.round(sig.reshape(K, N) * w), -32768, 32767).astype(np.int64)
    ref = np.mean(np.abs(np.fft.fft(q, axis=1)[:, :F] / N) ** 2, axis=0)
    c = int(round(f0)); far = np.abs(np.arange(F) - c) >= 40; far[:4] = False
    br = np.array([int(format(u, f'0{L}b')[::-1], 2) for u in range(N)])
    out = []
    for name, kw, sc in (('16-bit', dict(pre=0, W=16), 1.0), ('25-bit Q14', dict(pre=8, W=None), 2.0 ** -8), ('25-bit Q15', dict(pre=8, W=None, tw=15), 2.0 ** -8)):
        P = np.zeros(F)
        for k in range(0, K, 2):
            zr, zi = model(q[k], q[k + 1], L, **kw)
            Z = np.zeros(N, complex); Z[br] = (zr + 1j * zi) * sc
            P += (np.abs(Z[:F]) ** 2 + np.abs(Z[(-np.arange(F)) % N]) ** 2) / 2
        P /= K
        out.append('%s: median %+.2f dB, mean %+.2f dB, worst far bin %+.2f dB above the reference worst' % (
            name, 10 * np.log10(np.median(P[far]) / np.median(ref[far])), 10 * np.log10(np.mean(P[far]) / np.mean(ref[far])),
            10 * np.log10(P[far].max() / ref[far].max())))
    print('N=%d K=%d floor %d dBFS | ' % (N, K, fl) + ' | '.join(out))
for L, K in ((13, 16),):
    for fl in (-90, -100):
        run(L, K, fl)
