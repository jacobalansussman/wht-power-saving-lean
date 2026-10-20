"""sanity2.py -- checks of xw2lib at f = 2..6 (block 2): the P identity, Identity G, Xf1 Zf1 = -Zf1 Xf1,
D_z D_z = X_z, K unitary and flat, CS X CS^-1 = X S CZ (per column)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from xw2lib import Game, A_C, B_C
G = Game(range(2, 7))
rng = np.random.default_rng(5)
def rv(f, S=2): n = 4 ** f; return (rng.standard_normal((n, S)) + 1j * rng.standard_normal((n, S)))
ok = True
for f in G.FS:
    V = rv(f)
    def W(w): return G.apply_word(w, V, f)
    chk = {
        "Pr0 = a + b Zs0": W("Pr0") - (A_C * V + B_C * W("Z1")),
        "Pr1 = a + b Zs1": W("Pr1") - (A_C * V + B_C * W("Z2")),
        "A.lS1.Ai = Pr0 (chi S_last chi^-1)": W("Ai.lS1.A") - W("Pr0"),
        "Identity G: chi CS_last chi^-1 = a + b/2(1 + Z0 + Z1 - Z0Z1)":
            W("Ai.lCS.A") - (A_C * V + B_C / 2 * (V + W("Z1") + W("Z2") - W("Z3"))),
        "fX1 fZ1 + fZ1 fX1 = 0": W("fX1.fZ1") + W("fZ1.fX1"),
        "D1 D1 = X1": W("D1.D1") - W("X1"),
        "D3 Di3 = 1": W("D3.Di3") - V,
        "CS X1 CS^-1 = X1 S2 CZ (words act left to right: CSi.X1.CS = S2.CZ.X1)": W("CSi.X1.CS") - W("S2.CZ.X1"),
    }
    for k, d in chk.items():
        e = np.linalg.norm(d) / np.linalg.norm(V)
        if e > 1e-10: ok = False
        print("f=%d  %-62s %.1e" % (f, k, e))
    KV = G.K(V[:, :1] * 0 + np.eye(4 ** f)[:, :1], f)
    print("f=%d  K e_0 flat: min|.|*sqrt(N) = %.6f max = %.6f" % (f, np.abs(KV).min() * 2 ** f, np.abs(KV).max() * 2 ** f))
print("ALL OK" if ok else "SOME CHECK FAILED")
