#!/usr/bin/env python3
"""dr.py -- av-5: the dynamic-range test of at-cross H4 (wave-14/at-cross/h3/bfp.py), run on THIS design's bit-true
model (gen24.model: 25-bit words, exact, no wrap) and on as-6's 16-bit halving datapath (same model with W=16, pre=0,
wrap), as packed pairs (two Hann-windowed real records per complex FFT, PSD = (|Z_k|^2 + |Z_{N-k}|^2)/2, averaged).
Reference: float FFT of the same integer windowed frames.  Also the hardware PSD path: |Z|^2 of the integer outputs
accumulated exactly (48-bit integer) -- identical to the float of the same integers (all sums < 2^53)."""
import sys, numpy as np
sys.path.insert(0, '.')
from gen24 import model
rng = np.random.default_rng(5)
def run(L, K, fl, tone_db=-6.0):
    N = 1 << L; F = N // 2 + 1; t = np.arange(N * K); w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(N) / N)
    sig = 32767 * 10 ** (tone_db / 20) * np.cos(2 * np.pi * (N / 12.3) * t / N) + rng.normal(0, 32767 * 10 ** (fl / 20), N * K)
    q = np.clip(np.round(sig.reshape(K, N) * w), -32768, 32767).astype(np.int64)
    ref = np.mean(np.abs(np.fft.fft(q, axis=1)[:, :F] / N) ** 2, axis=0)
    c = int(round(N / 12.3)); msk = np.ones(F, bool); msk[max(c - 8, 0):c + 9] = False; msk[:4] = False
    res = {}
    for name, kw, scale in (('16-bit halving (as-6)', dict(pre=0, W=16), 1.0), ('25-bit (this design)', dict(pre=8, W=None), 2.0 ** -8)):
        P = np.zeros(F); acc = np.zeros(F, dtype=object) if name.startswith('25') else None
        maxacc = 0
        for k in range(0, K, 2):
            zr, zi = model(q[k], q[k + 1], L, **kw)
            br = np.array([int(format(u, f'0{L}b')[::-1], 2) for u in range(N)])
            Zr = np.zeros(N, np.int64); Zi = np.zeros(N, np.int64); Zr[br] = zr; Zi[br] = zi
            p2 = Zr ** 2 + Zi ** 2                                   # |Z_k|^2, exact integers (< 2^49)
            pk = p2[:F] + p2[(-np.arange(F)) % N]                    # |Z_k|^2 + |Z_{N-k}|^2
            maxacc = max(maxacc, int(pk.max()))
            P += pk.astype(float) * (scale ** 2) / 2 * (N ** 0)       # in units of the halved output
        P /= K; P *= N ** 0
        # the FFT output is X/N for both (every stage halves): compare shapes, the absolute scale cancels in the ratio
        Pn = P / np.median(P[~msk][:1]) if False else P
        res[name] = 10 * np.log10(np.median(Pn[msk]) / np.median(ref[msk]))
        if name.startswith('25'):
            res['acc_bits_needed (K/2 sums of |Z_k|^2+|Z_N-k|^2)'] = int(maxacc * (K // 2)).bit_length()
    truth = 10 * np.log10(np.median(ref[msk]) / ref.max())
    return truth, res
for L, K in ((13, 16), (13, 32), (14, 16)):
    for fl in (-70, -90, -100):
        truth, res = run(L, K, fl)
        print(f'N={1<<L:6d} K={K:3d} tone -6 dBFS + noise {fl} dBFS (truth floor {truth:.1f} dB below the tone bin) | ' + ' | '.join(f'{k}: {v:+.2f} dB' if isinstance(v, float) else f'{k}: {v}' for k, v in res.items()))
