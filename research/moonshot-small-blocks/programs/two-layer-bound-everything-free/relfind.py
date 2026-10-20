"""relfind.py -- mw2-xwlayers.  Exact-in-spirit relation finder for the XW family algebra at block 2.

Words with 0, 1 or 2 paid letters (plus truncated orbit words D1 g^j D2 of exponential-order g) are evaluated on S
random complex vectors at f = 2..6 jointly (ONE coefficient vector for all f; no dense operators).  Identical group
elements are merged first (fingerprints).  The kernel of the stacked evaluation matrix is the space of every-f relations
among the menu (up to the numerical cut, reported with its singular-value gap).  Then:
  (a) paid-count homogeneity: dim N  vs  sum over paid counts d of dim N_d   (excess = 'collapses')
  (b) layer homogeneity:      dim N  vs  sum over layer bins of dim N_bin   (excess = inhomogeneous relations)
  layer of a word = slope of log2(support of w e_x) between f = 4 and f = 6 (max over 3 basis vectors), rounded.
Every relation outside the homogeneous span is printed (largest terms) and re-tested at f = 7, 8 with fresh vectors.
Usage: relfind.py [menu]   menu in {main, orbit}
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, itertools, time
from xw2lib import Game

FS = [2, 3, 4, 5, 6]
S = 2
kind = sys.argv[1] if len(sys.argv) > 1 else "main"
t0 = time.time()

FREE = ["1", "P", "C01", "C10", "X1", "X2", "X3", "Z1", "Z2", "Z3", "S1", "S2", "S3", "CZ", "CS", "CSi", "Neg",
        "A", "Ai", "a0", "a1", "a0i", "a1i", "E01", "E10", "B", "Bi", "Pr0", "Pr1",
        "lX1", "lX2", "lZ1", "lZ2", "lS1", "lS2", "lCZ", "lCS", "fX1", "fX2", "fZ1", "fZ2", "fS1", "fS2", "fCZ", "fCS",
        "Ai.lS1.A", "Ai.lS2.A", "Ai.lCS.A", "A.CS.Ai", "Ai.CS.A", "B.CS.Bi", "A.B", "A.CS.B", "B.A", "Ai.fX1.A",
        "lS1.lS2", "Pr0.Pr1", "fS1.lS1"]
LFR = ["1", "CS", "A", "Pr0", "lCS", "Ai", "B", "a0"]
RFR = ["1", "CS", "Ai", "Pr1", "fCS", "A", "Bi", "a0i"]
MID = ["1", "CS", "CSi", "A", "Ai", "Pr0", "Pr1", "lCS", "fCS", "lS1", "A.CS.Ai", "B", "Z1", "X2", "S1", "CZ", "Ai.lCS.A",
       "fZ1", "lX1", "E01"]
words = list(FREE)
words += ["%s.D%d.%s" % (L, z, R) for L in LFR for z in (1, 2, 3) for R in RFR]
words += ["D%d.%s.D%d" % (z, Mi, w) for Mi in MID for z in (1, 2, 3) for w in (1, 2, 3)]
words += ["%s.D%d.%s.D%d.%s" % (L, z, Mi, w, R) for (L, R) in (("A", "Ai"), ("Pr0", "Pr1"), ("CS", "CSi"))
          for Mi in ("1", "CS", "Pr0") for z, w in ((1, 2), (1, 1), (2, 3))]
if kind == "ctl":   # positive control for the inhomogeneity detector: the sandwich h P0 h^-1 = a + b h Zs0 h^-1
    h, hi = "D1.CS.D2", "Di2.CSi.Di1"          # matrix h = D2 CS D1 (words act left to right)
    words = ["1", "Z1", "Pr0", "Z2", "X1", hi + ".Z1." + h, hi + ".Pr0." + h, hi + ".Z2." + h, hi + ".Pr1." + h,
             "D1.D2", "D1.CS.D2", "D2.Pr0.D1", "D2.Z1.D1"]
if kind == "orbit":
    for g in ("A.B", "A.CS.B", "B.A"):
        for j in range(0, 9):
            words.append("D1." + ".".join([g] * j + ["1"]) + ".D2")
            words.append("D3." + ".".join([g] * j + ["1"]) + ".D3")
G = Game(FS + [7, 8])
rng = np.random.default_rng(2026)
V = {f: (rng.standard_normal((4 ** f, S)) + 1j * rng.standard_normal((4 ** f, S))) / np.sqrt(2) for f in FS + [7, 8]}
def paid(w): return sum(1 for t in w.split(".") if t[:1] == "D" and t.lstrip("Di").isdigit())
def evalw(w, fs): return np.concatenate([G.apply_word(w, V[f], f).reshape(-1) for f in fs])
cols = []; keep = []; seen = {}; dup = 0
for w in words:
    c = evalw(w, FS)
    key = tuple(np.round(c[:64] * 1e6).astype(np.int64).view(np.int64).tolist()) if False else None
    h = hash(np.round(c, 8).tobytes())
    if h in seen and np.allclose(cols[seen[h]], c, atol=1e-9):
        dup += 1; continue
    seen[h] = len(cols); cols.append(c); keep.append(w)
Mt = np.stack(cols, axis=1)
print("menu %s: %d words, %d distinct group elements (%d duplicates merged); rows %d; %.1f s"
      % (kind, len(words), len(keep), dup, Mt.shape[0], time.time() - t0), flush=True)
# layer of each word
def layer(w):
    ss = []
    for f in (4, 6):
        n = 4 ** f; best = 0
        for x in (0, 5, n - 3):
            e = np.zeros((n, 1), complex); e[x] = 1
            best = max(best, int((np.abs(G.apply_word(w, e, f)) > 1e-12).sum()))
        ss.append(np.log2(best))
    return (ss[1] - ss[0]) / 2
lay = np.array([layer(w) for w in keep]); pc = np.array([paid(w) for w in keep])
lbin = np.where(lay > 1.5, 2, np.where(lay > 0.5, 1, 0))
print("layer slopes seen (rounded to 0.05):", sorted(set(np.round(lay * 20) / 20)), flush=True)
def nullspace(A, rel=1e-9):
    if A.shape[1] == 0: return np.zeros((0, 0)), None
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    k = int((s > rel * s[0]).sum())
    gap = (s[k - 1], s[k] if k < len(s) else 0.0)
    return Vh[k:].conj().T, gap
Mn = Mt / np.linalg.norm(Mt, axis=0)
N, gap = nullspace(Mn)
print("dim of relation space N = %d   (cut between singular values %.2e and %.2e)" % (N.shape[1], gap[0], gap[1]), flush=True)
def homog(bins, name):
    tot = 0; parts = []
    for bv in sorted(set(bins)):
        idx = np.where(bins == bv)[0]
        Nb, _ = nullspace(Mn[:, idx]); tot += Nb.shape[1]
        if Nb.shape[1]:
            full = np.zeros((len(keep), Nb.shape[1]), complex); full[idx] = Nb; parts.append(full)
    print("%s: sum of dims within bins = %d  ->  excess (relations not split by %s) = %d"
          % (name, tot, name, N.shape[1] - tot), flush=True)
    return np.concatenate(parts, axis=1) if parts else np.zeros((len(keep), 0))
Hp = homog(pc, "paid count")
Hl = homog(lbin, "layer bin")
def report(H, name, maxshow=12):
    if N.shape[1] == 0: return
    if H.shape[1]:
        Uh, sh, _ = np.linalg.svd(H, full_matrices=False); Q = Uh[:, sh > 1e-8 * sh[0]]
    else:
        Q = np.zeros((len(keep), 0))
    R = N - Q @ (Q.conj().T @ N)
    U, s, Vh = np.linalg.svd(R, full_matrices=False)
    k = int((s > 1e-6).sum())
    print("--- %d independent relations not split by %s; representatives:" % (k, name))
    for t in range(min(k, maxshow)):
        c = U[:, t]
        # sparsify: least-norm combination is not sparse; show top terms
        order = np.argsort(-np.abs(c))
        sup = [i for i in order if abs(c[i]) > 1e-6 * abs(c[order[0]])]
        c = c / c[order[0]]
        ver = {}
        for f in (7, 8):
            tot = sum(c[i] * G.apply_word(keep[i], V[f], f) for i in sup)
            nrm = sum(abs(c[i]) * np.linalg.norm(V[f]) for i in sup)
            ver[f] = np.linalg.norm(tot) / nrm
        print("  relation %d: %d terms; recheck rel. residual f=7: %.1e  f=8: %.1e" % (t, len(sup), ver[7], ver[8]))
        for i in sup[:10]:
            print("     %+.4f%+.4fi  [paid %d, layer %.2f]  %s" % (c[i].real, c[i].imag, pc[i], lay[i], keep[i]))
    sys.stdout.flush()
report(Hp, "paid count")
report(Hl, "layer bin")
print("done %.1f s" % (time.time() - t0))
