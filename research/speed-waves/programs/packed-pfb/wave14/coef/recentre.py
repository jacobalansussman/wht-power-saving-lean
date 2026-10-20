# Re-centre a prototype so that it is a true even-length linear-phase filter, h[n] = h[L-1-n] (critic ar-3 found that
# ap-4's prototypes are symmetric about index 7,160 CIRCULARLY). The FFT-interpolated prototype is band-limited and
# periodic, so a circular shift by 7.5 samples (phase ramp) is exact: centre 7,160 -> 7,167.5 = (L-1)/2.
import sys, numpy as np
src, dst = sys.argv[1], sys.argv[2]
h = np.load(src); L = h.size
c0 = 7160.0; sh = (L - 1) / 2 - c0
k = np.fft.fftfreq(L) * L
X = np.fft.fft(h) * np.exp(-2j * np.pi * k * sh / L)
# keep the Nyquist bin real (L even): a half-sample shift there makes it imaginary; it is ~0 for this lowpass anyway
hs = np.real(np.fft.ifft(X)); hs /= np.abs(hs).max()
print("source circular symmetry error about 7160:", np.abs(h - np.roll(h[::-1], 2 * 7160 + 1 - L)).max())
print("new: max |h[n] - h[L-1-n]| =", np.abs(hs - hs[::-1]).max(), " imag part dropped:", np.abs(np.fft.ifft(X).imag).max())
np.save(dst, hs)
