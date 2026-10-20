#!/usr/bin/env python3
"""thresholds.py (rm2-phase): what a Theorem-R certificate (= a frame family, = a test in Lemma N's block) must reach to
refuse the family's shapes W in 2^m - 1, with phases.  Exact rationals.
  need(W)    = (2^m - 1)/W                 (the certificate's ratio abar_m/abar_1 must EXCEED this)
  T_m(m)     = 2 - 2^(1-m)                 (spectral ceiling at depth m; exact for m <= 40, a2-vi-stability)
  angle      = need / T_m(m) = 2/c          (c = W / 2^(m-2), the sharing fraction): the angle ratio s_1/s_m needed
  lift lam   = m/need - 1                   (a lift proves m/(1 + lam); it must have edge excess per edge below this)
Known: best proved bounds 2.8, 41/13, 834241/278850, 1188222/399685 at m = 3..6; 2.9447 from m = 4 on (engine).
Hooks: value = T_m(m) * (angle ratio <= 2n/(n+2) < 2), so a hook never exceeds 2 T_m(m) = (2^m - 1)/2^(m-2).
Lifts: excess per edge >= 0.075 from block 4 on (meta2, upheld), so a lift never exceeds m/1.075."""
from fractions import Fraction as Fr
from math import log
best = {3: Fr(14, 5), 4: Fr(41, 13), 5: Fr(834241, 278850), 6: Fr(1188222, 399685)}
for m in range(3, 11):
    M = 2 ** m - 1; full = 2 ** (m - 2); T = 2 - Fr(1, 2 ** (m - 1))
    bar = m ** 0.9
    print(f"block {m}: M = {M}, full sharing W = {full}; T_m(m) = {T}; 2 T_m(m) = {2*T} = M/full: {2*T == Fr(M, full)}; "
          f"best proved bound {float(best.get(m, Fr(29447, 10000))):.4f}; bar allows {bar:.3f} per array")
    lo = int(M / bar) + 1
    for W in range(min(max(lo, 1), full), full + 3):
        need = Fr(M, W); s = log(m * W / M) / log(m)
        if s <= 0: continue
        c = Fr(W, full); ang = need / T; lam = Fr(m) / need - 1
        refused = best.get(m, Fr(29447, 10000)) * W > M
        print(f"   {W:3d} in {M}: saving {s:.3f}; per array {float(need):.4f}; sharing {float(c):.3f}; angle ratio needed "
              f"> {float(ang):.4f}; lift excess per edge needed < {float(lam):.4f}"
              f"{'; lifts can never (floor 0.075)' if lam <= Fr(3, 40) and m >= 4 else ''}"
              f"{'; REFUSED by the proved bound' if refused else ''}")
