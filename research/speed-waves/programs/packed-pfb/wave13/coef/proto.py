# SKA-Low-like prototype filter for the oversampled PFB (N = 1024, M = 864, 14 taps per branch).
# Method as in Comoretto et al. (arXiv 2002.10307, sec. 4.3): an equiripple (Remez) filter of 1/16 the size
# (896 taps, 64 branches), then FFT interpolation by 16 to 14,336 taps.
# Targets (paper): pass-band ripple +-0.17 dB over |f| <= 0.5 channel; > 60 dB at the start of the aliased
# region (0.685 channel = 32/27 - 1/2) and ~86 dB "for most of the stop-band".
# Writes proto_h.npy (float64, 14,336 taps, peak 1) and prints its response.
import sys, numpy as np
from scipy.signal import remez, resample
N, T, R = 1024, 14, 16
L = N * T; L0 = L // R; N0 = N // R
fp, fs1, fs2 = 0.5, 0.685, 1.5           # channel units
w1 = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
w2 = float(sys.argv[2]) if len(sys.argv) > 2 else 400.0
ch0 = 1.0 / N0
h0 = remez(L0, [0, fp*ch0, fs1*ch0, fs2*ch0, fs2*ch0 + 1e-6, 0.5], [1, 0, 0], weight=[1, w1, w2], fs=1.0, maxiter=200, grid_density=32)
h = resample(h0, L)                      # periodic FFT interpolation by 16
h = h / np.abs(h).max()

def resp_db(x, NF=1 << 21):
    H = np.abs(np.fft.rfft(x, NF)); H /= H[0]
    f = np.arange(NF // 2 + 1) / NF * N  # channel units
    return f, 20 * np.log10(np.maximum(H, 1e-30))

def report(x, tag):
    f, R_ = resp_db(x)
    pb = R_[f <= fp]; s1 = R_[(f >= fs1) & (f < fs2)]; s2 = R_[f >= fs2]
    rip = (pb.max() - pb.min()) / 2
    print(f"{tag}: ripple +-{rip:.3f} dB  worst stop [0.685,1.5) ch {s1.max():7.2f} dB  worst stop >=1.5 ch {s2.max():7.2f} dB"
          f"  median >=1.5 {np.median(s2):7.1f}")
    return rip, s1.max(), s2.max()

if __name__ == '__main__':
    report(h0 / np.abs(h0).max(), f"short 896 taps w1={w1} w2={w2}")
    report(h, f"interp 14336")
    np.save(__file__.replace('proto.py', 'proto_h.npy'), h)
