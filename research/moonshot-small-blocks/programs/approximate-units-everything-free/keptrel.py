"""keptrel.py -- xw-border.  Which relations does a FIBREWISE sub-call keep?  (floating point, random vectors)
C0  sanity: the fibrewise model with Y = Q reproduces xw2lib's exact D1, D2, D3.
C1  [D1, D2] = 0 for INDEPENDENT random sub-call values Y, Y' (a relation of degree 2 in the paid letters).
C2  commutation of D1 (random Y) with row-1-only letters: a1, a1i (running XOR on row 1: unbounded order), X2, Z2,
    lX2, lZ2, Pr1; of D2 (random Y') with row-0-only letters a0, X1, Z1, lZ1, lX1, Pr0.
C3  relations that a GENERIC fibrewise sub-call breaks: [D1, X1] (X1 = X on row 0 of every column), the twist
    lZ1 D1 lZ1 = i lX1 D1, D1 D1 = X1, and [D1, D3].
C4  w3-border's exotic relation R7 (pi- u - i lZ1 lX1 u_b) with (a) exact Q, (b) generic random Y, (c) Y equivariant
    under X^(x)f and exact on the last bit.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from keptlib import Fib

rng = np.random.default_rng(1)
rel = lambda a, b: float(np.linalg.norm(a - b) / max(np.linalg.norm(a), np.linalg.norm(b), 1e-300))
for f in range(1, 6):
    F = Fib(f)
    v = rng.standard_normal((F.n, 3)) + 1j * rng.standard_normal((F.n, 3))
    out = ["f=%d" % f]
    c0 = max(rel(F.D(z, v), F.G.apply_D(z, v, f)) for z in (1, 2, 3))
    out.append("C0 %.1e" % c0)
    Y, Yp = F.rand_Y(rng), F.rand_Y(rng)
    out.append("C1 [D1(Y),D2(Y')] %.1e" % rel(F.D(1, F.D(2, v, Yp), Y), F.D(2, F.D(1, v, Y), Yp)))
    c2 = []
    for w in ["a1", "a1i", "X2", "Z2", "lX2", "lZ2", "Pr1"]:
        c2.append(rel(F.D(1, F.free(w, v), Y), F.free(w, F.D(1, v, Y))))
    for w in ["a0", "X1", "Z1", "lZ1", "lX1", "Pr0"]:
        c2.append(rel(F.D(2, F.free(w, v), Yp), F.free(w, F.D(2, v, Yp))))
    out.append("C2 max %.1e" % max(c2))
    if f >= 2:
        out.append("C2' a1 order-growth control [a1, D2(Y')] %.2f" % rel(F.D(2, F.free("a1", v), Yp), F.free("a1", F.D(2, v, Yp))))
    c3 = {"[D1,X1]": rel(F.D(1, F.free("X1", v), Y), F.free("X1", F.D(1, v, Y))),
          "twist": rel(F.free("lZ1", F.D(1, F.free("lZ1", v), Y)), 1j * F.free("lX1", F.D(1, v, Y))),
          "D1D1=X1": rel(F.D(1, F.D(1, v, Y), Y), F.free("X1", v)),
          "[D1,D3]": rel(F.D(1, F.D(3, v, Yp), Y), F.D(3, F.D(1, v, Y), Yp))}
    c3q = {"[D1,X1]": rel(F.D(1, F.free("X1", v)), F.free("X1", F.D(1, v))),
           "twist": rel(F.free("lZ1", F.D(1, F.free("lZ1", v))), 1j * F.free("lX1", F.D(1, v))),
           "D1D1=X1": rel(F.D(1, F.D(1, v)), F.free("X1", v)),
           "[D1,D3]": rel(F.D(1, F.D(3, v)), F.D(3, F.D(1, v)))}
    out.append("C3 generic " + " ".join("%s %.2f" % kv for kv in c3.items()))
    out.append("C3 exact-Q " + " ".join("%s %.0e" % kv for kv in c3q.items()))
    # R7
    def R7(Y1, Y2):
        t = 0.5 * (v + F.free("X1", v))
        ua = F.D(2, t, Y2)
        u = F.D(1, t + F.free("lZ1", ua), Y1)
        ub = F.D(2, 0.5 * (u + F.free("X1", u)), Y2)
        pim = 0.5 * (u - F.free("X1", u))
        # word "lZ1.lX1" in xw2lib acts lZ1 first; we need  i * lZ1 lX1 u_b  (lX1 first, then lZ1)
        sec = 1j * F.free("lZ1", F.free("lX1", ub))
        return rel(pim, sec), float(np.linalg.norm(pim) / np.linalg.norm(v))
    Ye = F.eq_last_exact_Y(rng)
    r7 = [R7(None, None), R7(Y, Yp), R7(Ye, Yp)]
    out.append("C4 R7 rel.diff exact %.1e | generic %.2f | equivariant+last-exact %.1e  (|pi- u|/|v| %.2f %.2f %.2f)"
               % (r7[0][0], r7[1][0], r7[2][0], r7[0][1], r7[1][1], r7[2][1]))
    print("\n   ".join(out), flush=True)
