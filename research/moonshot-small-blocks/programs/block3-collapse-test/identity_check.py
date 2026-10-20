"""xw3-extract: the explicit (b) identities behind the repeated-direction flags, checked at f = 3..6 on fresh inputs.
(I1) "1;2;4;4" (B2): with y = D2 D1 v,  a * D4 Pr0 D4 y  -  b * Z1 D4 Pr0 D4 y  =  i * X4 y   (a = (1+i)/2, b = (1-i)/2):
     D4 Pr0 D4 = Pr0 X4 (Pr0 commutes with D4, D4 D4 = X4) and Pr0 = a + b Z1, a^2 - b^2 = i.
(I2) K9 runs (C2/C4/C5, mix P01 between moves 1 and 2): D2 P01 D1 = P01 D1 D1 = P01 X1, so with M the top move and
     u = P01 X1 v:  a * M Pr0 u  -  b * Z1 M Pr0 u ... only when Z1 commutes with M; for M = D4: = i * M u (one move);
     for M = D7 (C4) Z1 does not commute with D7 and the identity is instead  a * M Pr0 u - b * M Z1 Pr0 u = i * M u."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xw3lib import Game, A_C as a, B_C as b
for f in (3, 4, 5, 6):
    G = Game([f], W=1, S=1, seed=900 + f); v = G.inp(0)
    y = G.D(2, G.D(1, v))
    lhs = {f: a * G.D(4, G.W_("Pr0", G.D(4, y)))[f] - b * G.W_("Z1", G.D(4, G.W_("Pr0", G.D(4, y))))[f]}
    rhs = 1j * G.W_("X4", y)[f]
    print(f"f={f} (I1): |lhs - rhs|/|rhs| = {np.linalg.norm(lhs[f] - rhs) / np.linalg.norm(rhs):.1e}")
    u0 = G.D(2, G.W_("P01", G.D(1, v))); u = G.W_("P01", G.W_("X1", v))   # P01 X1 v (X1 acts first)
    print(f"f={f} D2 P01 D1 v = P01 X1 v: {np.linalg.norm(u0[f] - u[f]) / np.linalg.norm(u[f]):.1e}")
    for z in (4, 7):
        l2 = a * G.D(z, G.W_("Pr0", u))[f] - b * G.D(z, G.W_("Z1", G.W_("Pr0", u)))[f]
        r2 = 1j * G.D(z, u)[f]
        print(f"f={f} (I2) top D{z}: |a D Pr0 u - b D Z1 Pr0 u - i D u| / |D u| = {np.linalg.norm(l2 - r2) / np.linalg.norm(r2):.1e}")
# (I3) the exact form of the C4 survivors (coefficient ratio Z1/Pr0 = -b, read off out_C4s601.txt): for readout g in
#      {1, P01, X2}:  g D7 Pr0 D2 P01 D1 v  -  b * g D7 Z1 D2 P01 D1 v  =  a * g D7 P01 X1 v  (one paid move; g D7 P01 X1 v
#      = (g P01) D7 X1 v is a lower term because P01 commutes with D7).
for f in (3, 4, 5, 6):
    G = Game([f], W=1, S=1, seed=950 + f); v = G.inp(0)
    w = G.D(2, G.W_("P01", G.D(1, v)))
    for g in ("1", "P01", "X2"):
        lhs = G.W_(g, G.D(7, G.W_("Pr0", w)))[f] - b * G.W_(g, G.D(7, G.W_("Z1", w)))[f]
        rhs = a * G.W_(g, G.D(7, G.W_("P01", G.W_("X1", v))))[f]
        rhs2 = a * G.W_(g, G.W_("P01", G.D(7, G.W_("X1", v))))[f]
        print(f"f={f} (I3) g={g}: {np.linalg.norm(lhs - rhs) / np.linalg.norm(rhs):.1e}; P01 commutes with D7: "
              f"{np.linalg.norm(rhs - rhs2) / np.linalg.norm(rhs):.1e}")
