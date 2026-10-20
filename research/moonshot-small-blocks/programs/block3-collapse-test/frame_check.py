"""xw3-extract: is the frame A.T2.Ai of run C5 ("1;2;4,A.T2.Ai") trivial on D4 (C D4 C^-1 = D4)?  And B / A frames of C1/C3."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xw3lib import Game
for f in (3, 4, 5, 6):
    G = Game([f], W=1, S=2, seed=11); v = G.inp(0)
    for (z, C) in ((4, "A.T2.Ai"), (1, "B"), (4, "A"), (4, "T2"), (4, "CS01")):
        d = np.linalg.norm(G.move(z, C, v)[f] - G.D(z, v)[f]) / np.linalg.norm(v[f])
        print(f"f={f}: |C D{z} C^-1 v - D{z} v| / |v| for C={C}: {d:.1e}")
