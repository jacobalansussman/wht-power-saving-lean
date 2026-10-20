"""tnorm.py -- for a tri_als solution: residual, and the size of the unconstrained cross-term T (out1's y-part),
relative to |K y|, at fresh probes."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from tri_als import Prog, MENUS
npz, pat, menu, target = sys.argv[1], [int(t) for t in sys.argv[2].split(",")], MENUS[sys.argv[3]], sys.argv[4]
FS = [int(t) for t in sys.argv[5].split(",")] if len(sys.argv) > 5 else [2, 3, 4]
d = np.load(npz, allow_pickle=True)
pr = Prog(FS, pat, menu, 6, 4242, target)
P = {k: d["|".join(map(str, k))] for k in pr.keys}
for f in FS:
    outs, _, _ = pr.forward(P, f)
    S = outs[0].shape[1] // 2
    T = outs[0][:, S:]; KY = pr.T[f][1][:, S:]
    print("f=%d  residual %.4f   |T y|/|K y| = %.3f   max|coef| %.2f" % (f, pr.score(P, [f]), np.linalg.norm(T) / np.linalg.norm(KY),
          max(np.abs(v).max() for v in P.values())))
