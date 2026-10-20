# c23_open_mixed.py -- ms-mixed: mixed block-3 units W1 kernels + W2 two-row targets in M moves that stay OPEN after
# the termwise-tightness theorem (which excludes every unit exactly tight on the hook: 2.8 W1 + 1.8 W2 = M).
# Saving = 1 - tau with W1 3^tau + W2 2^tau = M (tau-transfer).  Certificate fronts (ratio kernel, ratio two-row):
#   widened: hook (2.8, 1.8); L2xL3 J={1} (2.5152, 1.8485); L3xL4 J={0,2} (2.2766, 1.8936)   [literature, computed]
#   engine : hook (2.8, 1.8); 105 (2.45, 1.9)                                               [Theorem M, unchecked]
import math
fronts = {'widened': [(2.8, 1.8), (1424/752*0 + 2.5152, 1.8485), (2.2766, 89/47)],
          'engine': [(2.8, 1.8), (49/20, 38/20)]}
def tau(W1, W2, M):
    lo, hi = 0.0, 1.0
    for _ in range(80):
        t = (lo + hi) / 2
        if W1 * 3**t + W2 * 2**t > M: hi = t
        else: lo = t
    return lo
for game, fr in fronts.items():
    rows = []
    for W1 in range(0, 15):
        for W2 in range(0, 60):
            if W1 + W2 == 0: continue
            for M in range(1, 3 * W1 + 2 * W2):
                need = max(a * W1 + b * W2 for a, b in fr)
                if M + 1e-9 < need: continue                       # excluded by a certificate
                if abs(2.8 * W1 + 1.8 * W2 - M) < 1e-9 and W2 > 0: continue   # hook-tight with a two-row target: excluded (termwise theorem)
                s = 1 - tau(W1, W2, M)
                if s >= 0.1: rows.append((round(s, 4), W1, W2, M, round(M - 2.8 * W1 - 1.8 * W2, 3)))
                break                                              # smallest open M for this (W1, W2)
    rows.sort(reverse=True)
    print(f"{game}: open mixed units with saving >= 0.1 (saving, W1, W2, M, hook slack in moves), best first; smallest W first among the top:")
    for r in rows[:12]: print("   ", r)
    small = sorted([r for r in rows if r[1] + r[2] <= 14], key=lambda r: (r[1] + r[2], -r[0]))
    print("   with <= 14 arrays:", small[:8])
