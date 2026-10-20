"""orth.py -- mw2-xwlayers.  Numerical illustration of the analytic heart of PROOF.md (Lemma O):
for free words c, c' and any word w with at most m - 1 = 1 paid letters,  |tr((c K c')^dagger w)| / N <= C 2^(-f/2).
Words with 2 paid letters may or may not be orthogonal (top layer vs collapsed).  Exact traces at f = 2..5
(identity blocks), Hutchinson estimates at f = 6, 7 (128 complex Gaussian vectors; std about 1/sqrt(128 N) <= 1.4e-3)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, itertools
from xw2lib import Game
FS = [2, 3, 4, 5, 6, 7]
G = Game(FS)
CC = ["1", "CS", "A", "Pr0", "lCS", "Ai.CS.A", "P.CS"]                 # frames c, c' around K
W1 = ["D1", "CS.D1", "D1.CS", "A.D1.Ai", "Ai.D2.A", "Pr0.D3.Pr1", "lCS.D1.fCS", "B.D2.CS.Bi", "CS.D3.A.CS",
      "Pr0", "1", "A.CS", "lX1.D2.lS1", "E01.D1.E10"]                    # <= 1 paid letter
W2 = ["D1.D2", "D1.CS.D2", "D2.A.D1", "D1.D1", "D1.CS.D1", "D1.Pr0.D1", "D3.CS.D3", "A.D1.Ai.D2", "D1.Pr1.D2.CS"]
rng = np.random.default_rng(11)
def inner(q, w, f):
    """|tr(Q^dagger W)| / N for Q = word q, W = word w (both as matrices at f)."""
    n = 4 ** f
    if f <= 5:
        tot = 0j
        step = 1024
        for s in range(0, n, step):
            E = np.zeros((n, min(step, n - s)), complex); E[np.arange(s, s + E.shape[1]), np.arange(E.shape[1])] = 1
            tot += np.vdot(G.apply_word(q, E, f), G.apply_word(w, E, f))
        return abs(tot) / n
    S = 128
    V = (rng.standard_normal((n, S)) + 1j * rng.standard_normal((n, S))) / np.sqrt(2)
    return abs(np.vdot(G.apply_word(q, V, f), G.apply_word(w, V, f))) / (n * S)
out = []
for group, WL in (("<=1 paid letter", W1), ("2 paid letters", W2)):
    print("== words w with %s: max over frames c, c' of |tr((c K c')^+ w)|/N, and that times 2^(f/2)" % group)
    for w in WL:
        vals = []
        for f in FS:
            best = 0
            for c, c2 in itertools.product(CC, repeat=2):
                q = "%s.D1.D2.%s" % (c2, c)            # matrix c K c'  (c' acts first)
                best = max(best, inner(q, w, f))
            vals.append(best)
        print("  %-16s " % w + " ".join("%.4f" % v for v in vals) + "   x2^(f/2): " +
              " ".join("%.3f" % (v * 2 ** (f / 2)) for v, f in zip(vals, FS)), flush=True)
