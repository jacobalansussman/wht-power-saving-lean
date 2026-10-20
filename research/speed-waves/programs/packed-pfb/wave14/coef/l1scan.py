# L1 structure of prototypes A and B: per branch, per tap index; predicted noise of grouping schemes.
import numpy as np, sys
N, T = 1024, 14
for P in ('A', 'B'):
    h = np.load('proto_h_w20_400.npy' if P == 'A' else 'proto_h_B.npy')
    H = h.reshape(T, N)            # H[t, k] = h[k + N t]
    S = h.sum(); dc_branch = H.sum(0)
    L1 = np.abs(H).sum(0)
    print(P, "sum h", S, "S/N", S / N, "branch sum min/max", dc_branch.min(), dc_branch.max())
    print("  L1 per branch max", L1.max(), "at k", L1.argmax(), "mean", L1.mean(), "max/mean dc", L1.max() / (S / N))
    pt = np.abs(H).max(1)          # max over branches, per tap index
    print("  max |h| per tap index:", np.round(pt, 4))
    print("  max L1 over branches of tap pairs (t, 13-t):", [round(float(np.abs(H[[t, 13 - t]]).sum(0).max()), 4) for t in range(7)])
    for s in (19, 17):
        K = (1 << (s - 1)) - 1
        g = K / (128 * L1.max())
        noise = np.sqrt(h.size / 12) / (g * S)
        print(f"  s={s}: whole-branch scale {g:.1f}, peak coef {g*np.abs(h).max():.0f}, rms quant noise rel DC {20*np.log10(noise):.1f} dB, est worst (+9.5) {20*np.log10(noise)+9.5:.1f}")
