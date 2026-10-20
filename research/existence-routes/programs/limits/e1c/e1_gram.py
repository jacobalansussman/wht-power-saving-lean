"""e1_gram.py -- ex2-limits experiment E1 (14 Oct 2026).
Is the f-ultralimit trace of the everything-free game XW faithful modulo exact relations on free (monomial) words?

For a menu of FREE XW words at block 2 (letters from xw2lib.py, copied unchanged from ../../moonshot/mw2-xwlayers/),
compute the EXACT Gram matrix N*G_f(i,j) = tr(w_i^* w_j) (a Gaussian integer: monomial words with Z4 phases) at each f
given on the command line, and save it.  Words are first merged up to an f-dependent global phase (projective
families), judged at f = 3, 4, 5 jointly.
Usage: python3 e1_gram.py <menu> <f> [<f> ...]      menu in {ctl, main}
"""
import sys, os, time, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, scipy.sparse as sp
from xw2lib import Game

menu = sys.argv[1]; FS = [int(a) for a in sys.argv[2:]]
t0 = time.time()
FREE = ["1", "P", "C01", "C10", "X1", "X2", "X3", "Z1", "Z2", "Z3", "S1", "S2", "S3", "CZ", "CS", "CSi",
        "A", "Ai", "a0", "a1", "a0i", "a1i", "E01", "E10", "B", "Bi", "Pr0", "Pr1",
        "lX1", "lX2", "lZ1", "lZ2", "lS1", "lS2", "lCZ", "lCS", "fX1", "fX2", "fZ1", "fZ2", "fS1", "fS2", "fCZ", "fCS",
        "Ai.lS1.A", "Ai.lS2.A", "Ai.lCS.A", "A.CS.Ai", "Ai.CS.A", "B.CS.Bi", "A.B", "A.CS.B", "B.A", "Ai.fX1.A",
        "lS1.lS2", "Pr0.Pr1", "fS1.lS1"]                       # = relfind.py's FREE list without Neg
GEN = ["P", "C01", "X1", "Z1", "Z2", "S1", "S3", "CS", "CZ", "A", "Ai", "a0", "a1i", "E01", "B", "Bi",
       "Pr0", "Pr1", "lX1", "lZ1", "lS1", "lCS", "lC01", "fX1", "fZ2", "fCS"]
if menu == "ctl":      # column-uniform words only: the trace must be ((tr h_i^* h_j)/4)^f exactly
    words = ["1", "P", "C01", "X1", "Z1", "S1", "CS", "CZ", "S3", "P.CS", "C01.S1", "CS.X1", "X1.CS", "Z2.S3"]
elif menu == "main":
    words = list(FREE) + list(GEN) + ["%s.%s" % (a, b) for a in GEN for b in GEN]
elif menu == "chain3":  # triple products of chain / parity / boundary letters only
    G3 = ["A", "Ai", "a0", "E01", "B", "Pr0", "Pr1", "lX1", "lS1", "lCS", "fX1", "fCS", "CS"]
    words = ["1"] + G3 + ["%s.%s" % (a, b) for a in G3 for b in G3] + \
            ["%s.%s.%s" % (a, b, c) for a in G3 for b in G3 for c in G3]
else:
    raise SystemExit("menu?")

# projective dedupe at f = 3, 4, 5
DF = [3, 4, 5]
Gd = Game(DF)
keep = []; seen = {}
for w in words:
    ms = Gd.mono(w)
    key = []
    for f in DF:
        g = ms[f]
        dph = (g.ph - g.ph[0]) % 4                       # phase up to a global constant
        key.append(hash(g.sig.tobytes()) ^ hash(dph.astype(np.int8).tobytes()))
    key = tuple(key)
    if key in seen:
        continue
    seen[key] = w; keep.append(w)
n = len(keep)
print("menu %s: %d words -> %d projectively distinct families (judged at f = 3,4,5); %.1f s"
      % (menu, len(words), n, time.time() - t0), flush=True)
with open("e1_words_%s.txt" % menu, "w") as fh:
    fh.write("\n".join(keep) + "\n")
PH = np.array([1, 1j, -1, -1j])
G = Game(FS)
for f in FS:
    t1 = time.time()
    N = 4 ** f
    SIG = np.empty((n, N), np.int64); PHS = np.empty((n, N), np.int8)
    for i, w in enumerate(keep):
        g = G.mono(w)[f]
        # the word sends basis vector e_v to i^ph[v] e_sig[v]
        SIG[i] = g.sig; PHS[i] = g.ph
    G.cache.clear()
    v = np.arange(N, dtype=np.int64)
    comb = (v[None, :] * N + SIG).reshape(-1)
    uniq, inv = np.unique(comb, return_inverse=True)
    del comb, uniq
    rows = np.repeat(np.arange(n), N)
    E = sp.csr_matrix((PH[PHS.reshape(-1)], (rows, inv.reshape(-1))), shape=(n, int(inv.max()) + 1))
    del rows, inv, SIG, PHS
    GN = (E.conj() @ E.T).toarray()                    # N * G_f, exact Gaussian integers in floats
    np.save("e1_gram_%s_f%d.npy" % (menu, f), GN)
    gf = GN / N
    ev = np.linalg.eigvalsh((gf + gf.conj().T) / 2)
    nz = ev[ev > 1e-9]
    print("f=%d N=%d: rank %d of %d; exact null %d; smallest nonzero eig %.4e; #eig in (1e-9,1e-2]: %d; %.1f s"
          % (f, N, nz.size, n, n - nz.size, nz.min() if nz.size else 0, ((ev > 1e-9) & (ev <= 1e-2)).sum(),
             time.time() - t1), flush=True)
print("done %.1f s" % (time.time() - t0))
