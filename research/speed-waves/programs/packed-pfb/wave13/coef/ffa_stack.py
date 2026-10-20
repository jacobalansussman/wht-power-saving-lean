# Does ao-h1's fast-FIR (FFA) step across frames stack with the two-signals-per-DSP packing? [derived numerically]
# Per branch k the 14 taps h_k[t] = h[k + N t]; FFA-2 uses sub-filters H0 = even taps, H1 = odd taps (8-bit inputs),
# and H0 + H1 (7 taps) on the 9-bit pre-sums x0 + x1. Packing bounds (low field of s bits, offset trick):
#   8-bit inputs, s = 19:  max_k L1 * 128 <= 2^18 ;  9-bit inputs need a 9-bit high field: s + 9 <= 27 -> s = 18:
#   max_k L1(H0+H1) * 256 <= 2^17.  The common coefficient scale g is set by the tightest sub-filter.
import numpy as np
h = np.load('proto_h_w20_400.npy'); N, T = 1024, 14
H = h.reshape(T, N).T                       # H[k, t]
L1_full = np.abs(H).sum(1).max()
h0, h1 = H[:, 0::2], H[:, 1::2]
L1_0, L1_1, L1_s = np.abs(h0).sum(1).max(), np.abs(h1).sum(1).max(), np.abs(h0 + h1).sum(1).max()
g_full = 2**18 / 128 / L1_full
g_ffa = min(2**18 / 128 / L1_0, 2**18 / 128 / L1_1, 2**17 / 256 / L1_s)
print(f"max L1 per branch: full 14 taps {L1_full:.3f}; H0 {L1_0:.3f}; H1 {L1_1:.3f}; H0+H1 {L1_s:.3f} (peak tap 1.0)")
print(f"largest coefficient scale: packed full cascade g = {g_full:.0f}; packed FFA (s=19 for H0, H1; s=18 for H0+H1) g = {g_ffa:.0f}"
      f" -> {20*np.log10(g_ffa/g_full):+.1f} dB of coefficient precision")
