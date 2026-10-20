# quick: SQNR / spur of the FFT model for different W16 constant precisions and guard bits
import numpy as np, fft4k as F
rng = np.random.default_rng(7); N = F.N; n = np.arange(N); fr = 8
xr = np.clip(np.round(rng.normal(0, 8192, (fr, N))), -32767, 32767).astype(np.int64); xi = np.clip(np.round(rng.normal(0, 8192, (fr, N))), -32767, 32767).astype(np.int64)
ph = rng.uniform(0, 6.28, (fr, 1))
tr = np.round(0.891 * 32767 * np.cos(2 * np.pi * 333 * n / N + ph)).astype(np.int64); ti = np.round(0.891 * 32767 * np.sin(2 * np.pi * 333 * n / N + ph)).astype(np.int64)
def run(tag):
    out = []
    for a, b in ((xr, xi), (tr, ti)):
        ref = F.reference(a, b); yr, yi = F.model(a, b, 'V7'); y = yr + 1j * yi
        out.append(10 * np.log10(np.mean(np.abs(ref) ** 2) / np.mean(np.abs(y - ref) ** 2)))
        if a is tr:
            Y = np.zeros_like(y); Y[:, F.BR12] = y; P = np.abs(Y) ** 2; P[:, 330:337] = 0
            out.append(10 * np.log10((P.max(axis=1) / (np.abs(Y) ** 2).max(axis=1)).max()))
    print(tag, 'noise SQNR %.2f dB, tone SQNR %.2f dB, tone worst spur %.1f dBc' % tuple(out))
full = dict(F.CSD)
for G in (4, 2):
    F.G = G; F.CSD = full; run('G=%d digits 4/4/5' % G)
F.G = 2
F.CSD = {'c1': [(1, 0), (-1, 4), (-1, 6), (1, 9)], 's1': [(1, 1), (-1, 3), (1, 7)], 'c2': full['c2']}; run('G=2 digits 4/3/5')
F.CSD = {'c1': [(1, 0), (-1, 4), (-1, 6)], 's1': [(1, 1), (-1, 3), (1, 7)], 'c2': [(1, 0), (-1, 2), (-1, 5), (-1, 7)]}; run('G=2 digits 3/3/4')
