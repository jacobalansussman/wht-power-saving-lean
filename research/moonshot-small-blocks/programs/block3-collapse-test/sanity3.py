"""sanity3.py -- ms-xw3: letter identities at block 3 (numerical, f = 1..4)."""
import numpy as np
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xw3lib import Game, flat, A_C, B_C
G = Game([1, 2, 3, 4], W=1, S=1, seed=3)
v = G.inp(0)
def same(X, Y): return max(np.abs(X[f] - Y[f]).max() for f in X)
def W(w, X=None): return G.W_(w, v if X is None else X)
# 1. P identity: Pr0 = a + b Zstr_0
print("Pr0 = a + b Z1 :", same(W("Pr0"), {f: A_C * v[f] + B_C * W("Z1")[f] for f in v}))
# 2. Tof X0 Tof^-1 = X0 CNOT(1->2)
print("T2.X1.T2 = X1.C12 :", same(W("T2.X1.T2"), W("X1.C12")))
# 3. CS X0 CS^-1 = X0 S_1 CZ  (CS01 . X1 . CS01i  vs  X1 . S2 . CZ01 in some order)
for cand in ["X1.S2.CZ01", "S2.CZ01.X1", "X1.CZ01.S2.S2.S2", "CZ01.S2.X1"]:
    print("CS01.X1.CS01i =", cand, same(W("CS01.X1.CS01i"), W(cand)))
# 4. A then Ai = 1 ; B then Bi = 1 ; a0 a0i
print("A.Ai = 1:", same(W("A.Ai"), v), " B.Bi:", same(W("B.Bi"), v), " a1.a1i:", same(W("a1.a1i"), v))
# 5. chi-conjugated last-column phase = parity phase:  Ai . lS1 . A  (Ai first)  vs Pr0
print("A.lS1.Ai vs Pr0:", same(W("A.lS1.Ai"), W("Pr0")), "  Ai.lS1.A vs Pr0:", same(W("Ai.lS1.A"), W("Pr0")),
      " A.fS1.Ai:", same(W("A.fS1.Ai"), W("Pr0")), " Ai.fS1.A:", same(W("Ai.fS1.A"), W("Pr0")))
# 6. D_z^2 = X^z (column-uniform)
print("D1 D1 = X1:", same(G.D(1, G.D(1, v)), W("X1")))
# 7. K = D1 D2 D4 commute
print("D1D2D4 = D4D2D1:", same(G.D(1, G.D(2, G.D(4, v))), G.D(4, G.D(2, G.D(1, v)))))
# 8. move with frame: inverse consistency
print("move(1,'T2.A') then inverse frame:", same(G.move(1, "T2.A", G.move(1, "T2.A", v)), G.W_("T2.A", G.W_("X1", G.mono(G.winv("T2.A"), v)))))
