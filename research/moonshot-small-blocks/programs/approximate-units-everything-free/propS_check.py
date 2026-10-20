"""propS_check.py -- xw-border.  Numerical check of Proposition S: an R8-seeded border family (pattern 1,2,1; eps^-1
readouts exactly R8's; contents and linking coefficients c + eps*delta with random free-word combinations delta; random
eps^0 readouts) has, at exact Q, the eps -> 0 limit equal to the output of the exact 3-move program u1', u2', u3' of
PROOF.md section 6.  Also: the same family evaluated with GENERIC fibrewise sub-calls converges (kept-formal), and with
generic FULL paid matrices it blows up like 1/eps (relation-using in the strict sense)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from keptlib import Fib

WORDS = ["1", "X1", "X2", "Z1", "Z2", "a0", "a1", "CS", "C01", "P", "lZ1", "lX2"]
rng = np.random.default_rng(3)
for f in (2, 3, 4):
    F = Fib(f)
    x = rng.standard_normal((F.n, 2)) + 1j * rng.standard_normal((F.n, 2))
    def comb():
        c = rng.standard_normal(len(WORDS)) + 1j * rng.standard_normal(len(WORDS))
        return lambda v, c=c: sum(ci * F.free(w, v) for ci, w in zip(c, WORDS)) / len(WORDS)
    dc1, dc2, dc3, dB21, dB32, dB31, A1, A2, A3 = [comb() for _ in range(9)]
    pp = lambda v: 0.5 * (v + F.free("X1", v)); pm = lambda v: 0.5 * (v - F.free("X1", v))
    def border(eps, paid):
        u1 = paid(1, pp(x) + eps * dc1(x))
        u2 = paid(2, pp(x) + eps * dc2(x) + pm(u1) + eps * dB21(u1))
        u3 = paid(1, pp(u2) + eps * dB32(u2) + eps * dc3(x) + eps * dB31(u1))
        return A1(u1) + A2(u2) + A3(u3) + (pm(u2) - pm(u3)) / eps
    paidQ = lambda z, v: F.D(z, v)
    Y, Yp = F.rand_Y(rng), F.rand_Y(rng)
    paidY = lambda z, v: F.D(z, v, Y if z == 1 else Yp)
    G = {z: (rng.standard_normal((F.n, F.n)) + 1j * rng.standard_normal((F.n, F.n))) / np.sqrt(2 * F.n) for z in (1, 2)}
    paidG = lambda z, v: G[z] @ v
    # exact twin at Q
    u1p = F.D(1, pp(x) + pm(dc1(x)))
    u2p = F.D(2, pp(x) + pm(dc2(x)) + pm(u1p) + pm(dB21(pp(u1p))))
    u3p = F.D(1, pp(u2p) + pm(dB32(pp(u2p))) + pm(dc3(x)) + pm(dB31(pp(u1p))))
    twin = A1(pp(u1p)) + A2(pp(u2p)) + A3(pp(u3p)) + pm(u2p) - pm(u3p)
    nx = np.linalg.norm(x)
    line = ["f=%d" % f]
    for eps in (1e-2, 1e-4, 1e-6):
        bQ = border(eps, paidQ); bY = border(eps, paidY); bG = border(eps, paidG)
        line.append("eps=%.0e: |border_Q - twin|/|x| %.1e  |border_Y| %.2e  |border_G| %.2e" %
                    (eps, np.linalg.norm(bQ - twin) / nx, np.linalg.norm(bY) / nx, np.linalg.norm(bG) / nx))
    print("\n   ".join(line), flush=True)
