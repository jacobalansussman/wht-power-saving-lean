"""example_m3.py -- the simplest extra route, block size 3, checked exactly:
K = D(e3) (a + b t2)(a + b t1), t1 = direction e1 on the half x3 = 0 and direction e2 on the half x3 = 1
(t1 = F X(e1) F^-1 with F the Fredkin gate: swap bits 1,2 when bit 3 is set), t2 = the other way round.
Its first frame B(a + b t1) = B D(e1) F^-1 is the affine frame B D(e1) relabelled on the input side by F."""
import sys
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from rt import Game, A_, B_
g = Game(3); n = 8
def perm(f):
    M = np.zeros((n, n));
    for x in range(n): M[f(x), x] = 1
    return M
t1 = perm(lambda x: x ^ (1 if not (x & 4) else 2))
t2 = perm(lambda x: x ^ (2 if not (x & 4) else 1))
F = perm(lambda x: x if not (x & 4) else ((x & 4) | ((x & 1) << 1) | ((x & 2) >> 1)))
I = np.eye(n)
P = (A_ * I + B_ * g.Xm(4)) @ (A_ * I + B_ * t2) @ (A_ * I + B_ * t1)
print("K = D(e3)(a+b t2)(a+b t1) exactly:", np.allclose(P, g.K))
print("t1 = F X(e1) F^-1:", np.allclose(t1, F @ g.Xm(1) @ F.T))
print("t1, t2 commute:", np.allclose(t1 @ t2, t2 @ t1), "; t1 commutes with X(e3):", np.allclose(t1 @ g.Xm(4), g.Xm(4) @ t1))
print("frame B(a+b t1) == frame B D(e1) F^-1:", g.fkey(A_ * I + B_ * t1) == g.fkey(g.Dm(1) @ F.T))
print("frame of D(e2)D(e1) fixed by F^-1:", g.fkey(g.Dm(2) @ g.Dm(1) @ F.T) == g.fkey(g.Dm(2) @ g.Dm(1)))
print("F affine? (F(x)+F(y)+F(0) == F(x+y) for all x,y):", all((F.argmax(axis=0)[x] ^ F.argmax(axis=0)[y] ^ F.argmax(axis=0)[0]) == F.argmax(axis=0)[x ^ y] for x in range(n) for y in range(n)))
