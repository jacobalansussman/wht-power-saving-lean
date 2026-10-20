"""r8check.py -- xw-border.  A complexity-3 KEPT relation of degree 2 at block 2 that needs NO equivariance of the
sub-calls (unlike w3-border's R7): packed commutation through row-0 projections, which commute with every fibrewise
D2 sub-call.   pi+- = (1 +- X1)/2  (X1 = X on row 0 of every column; commutes with D2 = I (x) Y' for every Y').
  R8  (pattern 1,2,1):  u1 = D1 pi+ x,  u2 = D2(pi+ x + pi- u1),  u3 = D1(pi+ u2);   Psi = pi- u2 - pi- u3
        = pi- (D2 D1 - D1 D2) pi+ x   (uses only [D2, pi] = 0 and [D1, D2] = 0).
  R8' (mirror, pattern 2,1,2) with rho+- = (1 +- X2)/2 on row 1.
Checks at f = 1..5: Psi = 0 for exact Q and for random fibrewise Y, Y' (independent); Psi != 0 for random FULL
matrices in place of the paid letters (so it is nonzero in the strict formal algebra), and each half nonzero."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from keptlib import Fib

rng = np.random.default_rng(7)
for f in range(1, 6):
    F = Fib(f)
    x = rng.standard_normal((F.n, 2)) + 1j * rng.standard_normal((F.n, 2))
    full = {z: (rng.standard_normal((F.n, F.n)) + 1j * rng.standard_normal((F.n, F.n))) / np.sqrt(2 * F.n) for z in (1, 2)}
    Y, Yp = F.rand_Y(rng), F.rand_Y(rng)
    def run(paid, pz, qz, P):
        pp = lambda v: 0.5 * (v + F.free(P, v)); pm = lambda v: 0.5 * (v - F.free(P, v))
        u1 = paid(pz, pp(x))
        u2 = paid(qz, pp(x) + pm(u1))
        u3 = paid(pz, pp(u2))
        a, b = pm(u2), pm(u3)
        return np.linalg.norm(a - b) / np.linalg.norm(x), np.linalg.norm(a) / np.linalg.norm(x)
    res = []
    for (pz, qz, P) in [(1, 2, "X1"), (2, 1, "X2")]:
        ex = run(lambda z, v: F.D(z, v), pz, qz, P)
        fb = run(lambda z, v: F.D(z, v, Y if z == 1 else Yp), pz, qz, P)
        fu = run(lambda z, v: full[z] @ v, pz, qz, P)
        res.append("pattern %d,%d,%d: exact %.1e | fibrewise %.1e | full %.2f (half %.2f)" % (pz, qz, pz, ex[0], fb[0], fu[0], fu[1]))
    print("f=%d  " % f + "   ".join(res), flush=True)
