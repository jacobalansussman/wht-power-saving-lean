#!/usr/bin/env python3
"""H-A (an-cross, wave 12): am-1's precision lessons (floor-halving bias; Parseval-informed alternate halving)
applied to am-5's 16-bit streaming WHT arithmetic on the iCE40-UP5K.

Variants (all 16-bit words, radix-2 stages in SDF order, swapped butterfly irrelevant for error):
  V0  am-5 as built: every stage halves its inputs first (floor), then a-b / a+b  (scale 1/N)
  V1  every stage halves, but after the add, floor:   (a+-b)>>1                     (scale 1/N, 17-bit adder)
  V2  every stage halves after the add, round half up: (a+-b+1)>>1 (carry-in = 1)  (scale 1/N, 17-bit adder)
  V3  am-1's plan: only every other stage halves (after the add, round half up); the others add with
      saturation to 16 bits (scale 2^-ceil(L/2) ~ 1/sqrt(N)); counts saturated frames.
Reference: exact H x times the same scale, in float. Metric: SQNR (dB) over all outputs, max |err| in LSB,
mean error at output 0 (DC bias) in LSB. Inputs: frames of N 16-bit samples of several kinds.
"""
import sys
import numpy as np

def stages(v, L, mode):
    N = v.shape[1]
    sat = np.zeros(v.shape[0], dtype=bool)
    for i in range(L):
        D = N >> (i + 1)
        w = v.reshape(v.shape[0], -1, 2, D)
        a, b = w[:, :, 0, :], w[:, :, 1, :]
        if mode == 'V0':
            a, b = a >> 1, b >> 1
            s, d = a + b, a - b
        elif mode == 'V1':
            s, d = (a + b) >> 1, (a - b) >> 1
        elif mode == 'V2':
            s, d = (a + b + 1) >> 1, (a - b + 1) >> 1
        elif mode == 'V5':  # halve first (as am-5), constant carry-in 1 into the SUM on odd stages (SB_MAC16 BOTADDSUB_CARRYSELECT=1)
            a, b = a >> 1, b >> 1
            s, d = a + b + (i % 2), a - b
        elif mode == 'V6':  # ap-6: halve first, carry-in = the dropped LSB of the LATER input b (SB_MAC16 CI pin, BOTADDSUB_CARRYSELECT=3), every stage
            b0 = b & 1
            a, b = a >> 1, b >> 1
            s, d = a + b + b0, a - b
        elif mode == 'V4':
            if i % 2 == 0:
                s, d = (a + b) >> 1, (a - b) >> 1
            else:
                s, d = (a + b + 1) >> 1, (a - b + 1) >> 1
        elif mode == 'V3':
            if i % 2 == 0:
                s, d = (a + b + 1) >> 1, (a - b + 1) >> 1
            else:
                s, d = a + b, a - b
                over = (np.abs(s) > 32767).any(axis=(1, 2)) | (np.abs(d) > 32767).any(axis=(1, 2))
                sat |= over
                s, d = np.clip(s, -32768, 32767), np.clip(d, -32768, 32767)
        v = np.stack([s, d], axis=2).reshape(v.shape[0], N)
    return v, sat

def exact(x, L):
    v = x.astype(np.float64).copy()
    N = v.shape[1]
    for i in range(L):
        D = N >> (i + 1)
        w = v.reshape(v.shape[0], -1, 2, D)
        a, b = w[:, :, 0, :].copy(), w[:, :, 1, :].copy()
        v = np.stack([a + b, a - b], axis=2).reshape(v.shape[0], N)
    return v

def inputs(kind, F, N, rng):
    if kind == 'uniform16':
        x = rng.integers(-32768, 32768, size=(F, N))
    elif kind == 'gauss-12dBFS':
        x = np.round(rng.normal(0, 32768 / 4, size=(F, N)))
    elif kind == 'student3-20dBFS':
        x = np.round(rng.standard_t(3, size=(F, N)) * 32768 / 10)
    elif kind == 'tone+noise':
        n = np.arange(N)
        x = np.round(8000 * np.sin(2 * np.pi * 0.1234 * n + rng.uniform(0, 6.28, size=(F, 1))) + rng.normal(0, 800, size=(F, N)))
    elif kind == 'nonneg(image-like)':
        # single-pixel-imaging-like: non-negative intensities, smooth plus texture, mean ~ 1/4 full scale
        base = rng.uniform(4000, 12000, size=(F, 1))
        x = np.round(base + np.cumsum(rng.normal(0, 200, size=(F, N)), axis=1))
        x = np.clip(x, 0, 32767)
    return np.clip(x, -32768, 32767).astype(np.int64)

import os
MODES = os.environ.get('MODES', 'V0,V1,V2,V3').split(',')
HEAD = int(os.environ.get('HEAD', '0'))
LS = [int(t) for t in os.environ.get('LS', '10,12,14').split(',')]
def main():
    rng = np.random.default_rng(12)
    print('N kind variant SQNR_dB max|err|_LSB DCbias_LSB satFrames/F  (scale: V0-V2 1/N, V3 2^-ceil(L/2))')
    for L in LS:
        N = 1 << L
        F = max(8, 2 ** 18 // N)
        for kind in ('uniform16', 'gauss-12dBFS', 'student3-20dBFS', 'tone+noise', 'nonneg(image-like)'):
            x = inputs(kind, F, N, rng)
            H = exact(x, L)
            for mode in MODES:
                y, sat = stages(x.copy(), L, mode)
                halvings = L if mode != 'V3' else (L + 1) // 2
                ref = H / (1 << halvings)
                if mode == 'V3' and HEAD:
                    y, sat = stages(x.copy() >> HEAD, L, mode)
                    ref = H / (1 << (halvings + HEAD))
                ok = ~sat
                if ok.sum() == 0:
                    print(L and N, kind, mode, 'all frames saturate', '%d/%d' % (sat.sum(), F))
                    continue
                e = (y - ref)[ok]
                sq = 10 * np.log10((ref[ok] ** 2).mean() / (e ** 2).mean())
                print(N, kind, mode, '%.2f' % sq, '%.2f' % np.abs(e).max(), '%.2f' % (y[ok, 0] - ref[ok, 0]).mean(),
                      '%d/%d' % (sat.sum(), F))
            sys.stdout.flush()

if __name__ == '__main__':
    main()
