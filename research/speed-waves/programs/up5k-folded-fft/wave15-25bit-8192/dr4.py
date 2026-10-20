#!/usr/bin/env python3
"""dr4.py -- av-5: where the 25-bit design's worst far bin sits relative to the tone bin (spur-free range), -6 dBFS tone
plus -90 dBFS noise, 16 records of 8,192, packed pairs; and the worst far bin's position."""
import sys, numpy as np
sys.path.insert(0, '.')
from gen24 import model
rng = np.random.default_rng(11)
L, K, fl = 13, 16, -90
N = 1 << L; F = N // 2 + 1; t = np.arange(N * K); w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(N) / N); f0 = N / 12.3
sig = 32767 * 10 ** (-6 / 20) * np.cos(2 * np.pi * f0 * t / N) + rng.normal(0, 32767 * 10 ** (fl / 20), N * K)
q = np.clip(np.round(sig.reshape(K, N) * w), -32768, 32767).astype(np.int64)
ref = np.mean(np.abs(np.fft.fft(q, axis=1)[:, :F] / N) ** 2, axis=0)
c = int(round(f0)); far = np.abs(np.arange(F) - c) >= 40; far[:4] = False
br = np.array([int(format(u, f'0{L}b')[::-1], 2) for u in range(N)])
top = ref.max()
print('reference: worst far bin %.1f dB below the tone bin (bin %d); median floor %.1f dB below' % (
    10 * np.log10(top / ref[far].max()), np.arange(F)[far][ref[far].argmax()], 10 * np.log10(top / np.median(ref[far]))))
for name, kw, sc in (('16-bit', dict(pre=0, W=16), 1.0), ('25-bit Q14', dict(pre=8), 2.0 ** -8), ('25-bit Q15', dict(pre=8, tw=15), 2.0 ** -8)):
    P = np.zeros(F)
    for k in range(0, K, 2):
        zr, zi = model(q[k], q[k + 1], L, **kw)
        Z = np.zeros(N, complex); Z[br] = (zr + 1j * zi) * sc
        P += (np.abs(Z[:F]) ** 2 + np.abs(Z[(-np.arange(F)) % N]) ** 2) / 2
    P /= K
    fb = np.arange(F)[far]; o = np.argsort(P[far])[::-1][:5]
    print('%-11s worst far bins (dB below the tone bin @ bin): ' % name + ', '.join('%.1f @ %d' % (10 * np.log10(P.max() / P[far][i]), fb[i]) for i in o))
