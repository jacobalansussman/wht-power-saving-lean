"""v-xw3-extract: the (b) identities and the C5 frame, re-checked with my own letters (vlet.py), full rows, f = 3..7"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlet import P, A_ as a, B_ as b
def r(x, y): return np.linalg.norm(x - y) / np.linalg.norm(y)
for f in (3, 4, 5, 6, 7):
    Pp = P(f, 1, 4000 + f); v = Pp.v; W = Pp.word; D = Pp.D
    y = D(2, D(1, v))
    i1 = r(a * D(4, W("Pr0", D(4, y))) - b * W("Z1", D(4, W("Pr0", D(4, y)))), 1j * W("X4", y))
    dd = r(D(4, W("Pr0", D(4, y))), W("X4", W("Pr0", y)))            # D4 Pr0 D4 = Pr0 X4
    w = D(2, W("P01", D(1, v))); u = W("P01", W("X1", v))
    i2 = r(w, u)                                                       # D2 P01 D1 = P01 X1
    i3 = max(r(W(g, D(7, W("Pr0", w))) - b * W(g, D(7, W("Z1", w))), a * W(g, D(7, u))) for g in ("1", "P01", "X2"))
    # the RHS of I3 is a finite-menu lower word R9 . D7 . H30 . v for g in {1, P01, X2}
    low3 = max(r(W(g, D(7, u)), W(R, D(7, W(h, v)))) for g, R, h in (("1", "X2", "P01"), ("P01", "X1", "1"), ("X2", "1", "P01")))
    # C2 survivor 1: alpha (D4 Pr0 u) + i alpha (Z1 D4 Pr0 u) = alpha (a + i b) X2 D4 P01 v  (one move, menu word)
    c2 = r(D(4, W("Pr0", u)) + 1j * W("Z1", D(4, W("Pr0", u))), (a + 1j * b) * W("X2", D(4, W("P01", v))))
    # C5 frame: (A.T2.Ai) D4 (A.T2.Ai)^-1 = D4
    c5 = r(Pp.move(4, "A.T2.Ai", v), D(4, v))
    # small-f coincidence source: A^4 = 1 exactly when f <= 4
    A4 = r(W("A.A.A.A", v), v)
    print(f"f={f}: I1 {i1:.1e}; D4 Pr0 D4 = Pr0 X4 {dd:.1e}; D2 P01 D1 = P01 X1 {i2:.1e}; I3 {i3:.1e}; I3 rhs in menu {low3:.1e}; "
          f"C2 survivor-1 identity {c2:.1e}; C5 frame {c5:.1e}; |A^4 v - v|/|v| {A4:.2f}", flush=True)
