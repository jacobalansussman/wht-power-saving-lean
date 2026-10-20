"""rescore.py -- v-w7-xwborder.  Re-score saved w7-xwborder triangular solutions on NEW probe seeds (not 4242, not the
fit seeds) at f = 2..6, with w7-xwborder's own forward model (tri_als.Prog) -- a fresh-input check of the numbers in
summary.txt / RUNLOG, not an independent model.  Prints the unpenalised triangular residual and |T y| / |K y|.
usage: python3 -I rescore.py npz pattern menu target seed fcsv"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from tri_als import Prog, MENUS
npz, pat, menu, target, seed = sys.argv[1], [int(t) for t in sys.argv[2].split(",")], MENUS[sys.argv[3]], sys.argv[4], int(sys.argv[5])
FS = [int(t) for t in sys.argv[6].split(",")]
os.environ["MU"] = "0"
d = np.load(npz, allow_pickle=True)
pr = Prog(FS, pat, menu, 6, seed, target)
P = {k: d["|".join(map(str, k))] for k in pr.keys}
for f in FS:
    outs, _, _ = pr.forward(P, f)
    S = outs[0].shape[1] // 2
    T = outs[0][:, S:]; KY = pr.T[f][1][:, S:]
    tri = pr.score(P, [f])          # MU = 0 -> out1's y-part unconstrained: the triangular residual
    tn = np.linalg.norm(T) / np.linalg.norm(KY)
    print("%s f=%d  tri %.4f  |Ty|/|Ky| %.3f  product %.3f  max|coef| %.2f" % (os.path.basename(npz), f, tri, tn, tri * tn,
          max(np.abs(v).max() for v in P.values())), flush=True)
