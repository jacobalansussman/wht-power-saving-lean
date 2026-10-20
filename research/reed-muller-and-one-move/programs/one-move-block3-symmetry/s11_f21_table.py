#!/usr/bin/env python3
"""s11_f21_table.py (t2-block3, stdlib): F21-symmetric block-3 shapes (PROOF.md S9).

Array orbits: alpha of size 21, beta of size 7 (stabiliser Z3), gamma of size 3 (stabiliser Z7), delta fixed.
Move orbits: r of size 21, q of size 7 (stabiliser Z3; no move is fixed by an element of order 7, S3).
Screens:
  Z7 level (S4, S5):  80 j >= 224 a + 32 b,  j >= a + b + 2   with a = 3 alpha + beta, b = 3 gamma + delta, j = 3 r + q
  orbit count (S5'):  r + q >= alpha + beta + gamma + delta + 2
  Z3 characters of the Z7-trivial problem (S9): for k = 0, 1, 2
      80 r + q sigma_{-k} >= 224 alpha + d_k beta + 32 gamma + p_k delta,   d = (76,74,74), p = (12,10,10), sigma = (28,26,26)
      and r + q >= alpha + gamma + 2 for k != 0 when alpha + gamma >= 1 (frame-visible inputs only)
"""
import math
d, p, sg = (76, 74, 74), (12, 10, 10), (28, 26, 26)
lines = []
def survivors(W, M):
    out = []
    for al in range(W // 21 + 1):
        for be in range((W - 21 * al) // 7 + 1):
            for ga in range((W - 21 * al - 7 * be) // 3 + 1):
                de = W - 21 * al - 7 * be - 3 * ga
                a, b = 3 * al + be, 3 * ga + de
                for r in range(M // 21 + 1):
                    if (M - 21 * r) % 7:
                        continue
                    q = (M - 21 * r) // 7
                    j = 3 * r + q
                    if 80 * j < 224 * a + 32 * b or j < a + b + 2:
                        continue
                    if r + q < al + be + ga + de + 2:
                        continue
                    if any(80 * r + q * sg[(-k) % 3] < 224 * al + d[k] * be + 32 * ga + p[k] * de for k in range(3)):
                        continue
                    if al + ga >= 1 and r + q < al + ga + 2:
                        continue
                    out.append((al, be, ga, de, r, q, "=smaller unit" if b == 0 else "new"))
    return out
for W in range(5, 61):
    for M in range(math.ceil(2.8 * W - 1e-9), 3 * W):
        S = survivors(W, M)
        if S:
            s = "%d in %d (saving %.4f): %s" % (W, M, math.log(3 * W / M) / math.log(3),
                ", ".join("(arr 21x%d 7x%d 3x%d 1x%d; moves 21x%d 7x%d; %s)" % x for x in S[:3]) + (" ... [%d shapes]" % len(S) if len(S) > 3 else ""))
            print(s); lines.append(s)
if not lines:
    print("no F21-symmetric shape with W <= 60 saves a move"); lines.append("none")
open("out_s11_f21_table.txt", "w").write("\n".join(lines) + "\n")
