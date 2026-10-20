"""e1_paid.py -- ex2-limits E1b (14 Oct 2026).  Same question as e1_gram.py (is the XW f-limit trace faithful modulo
exact relations?) for words with ONE paid letter: w = L . D_z . R (letters act left to right, as in xw2lib.py: first
L, then D_z, then R; as a matrix w = R D_z L).  Gram G_f(i,j) = tr(w_i^* w_j)/N computed from sparse matrices
(D_z = dir(z)^(x)f has 2^f nonzeros per row).  Words are merged when equal up to a global phase at f = 3 and 4.
Usage: python3 e1_paid.py <f> [<f> ...]
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, scipy.sparse as sp
from xw2lib import Game, A_C, B_C, PH

FS = [int(a) for a in sys.argv[1:]]
t0 = time.time()
SIDE = ["1", "CS", "A", "Ai", "Pr0", "Pr1", "lCS", "fCS", "B", "a0", "S1", "C01", "lX1", "fX1", "E01", "P"]
words = ["%s.D%d.%s" % (L, z, R) for L in SIDE for z in (1, 2, 3) for R in SIDE]

def mono_mat(g):
    n = g.sig.size
    return sp.csr_matrix((PH[g.ph], (g.sig, np.arange(n))), shape=(n, n))

def dir_mat(z, f, m=2):
    q = 1 << m
    d = np.zeros((q, q), complex)
    for x in range(q):
        d[x, x] += A_C; d[x ^ z, x] += B_C
    M = sp.csr_matrix(d)
    out = sp.csr_matrix(np.ones((1, 1)))
    for j in range(f):                        # column j is bits m*j..m*j+m-1: kron with the new column on the left
        out = sp.kron(M, out, format="csr")
    return out

def word_mat(G, w, f):
    L, Dz, R = w.split(".D")[0], w.split(".D")[1].split(".")[0], ".".join(w.split(".D")[1].split(".")[1:])
    ML = mono_mat(G.mono(L)[f]); MR = mono_mat(G.mono(R)[f])
    return (MR @ dir_mat(int(Dz), f) @ ML).tocsr()

# dedupe up to a global phase, judged at f = 3, 4 on two random vectors
DF = [3, 4]
Gd = Game(DF)
rng = np.random.default_rng(7)
V = {f: rng.standard_normal((4 ** f, 2)) + 1j * rng.standard_normal((4 ** f, 2)) for f in DF}
for w in words[::97]:                         # control: sparse word matrices agree with xw2lib's apply_word
    for f in DF:
        err = np.abs(word_mat(Gd, w, f) @ V[f] - Gd.apply_word(w, V[f], f)).max()
        assert err < 1e-12, (w, f, err)
print("control: word matrices agree with xw2lib.apply_word on", len(words[::97]), "words at f = 3, 4", flush=True)
keep = []; reps = []
for w in words:
    fp = np.concatenate([(word_mat(Gd, w, f) @ V[f]).reshape(-1) for f in DF])
    dup = False
    for r in reps:
        c = np.vdot(r, fp) / np.vdot(r, r)
        if abs(abs(c) - 1) < 1e-9 and np.linalg.norm(fp - c * r) < 1e-8 * np.linalg.norm(fp):
            dup = True; break
    if not dup:
        reps.append(fp); keep.append(w)
n = len(keep)
print("paid menu: %d words -> %d distinct up to global phase (f = 3, 4); %.1f s" % (len(words), n, time.time() - t0),
      flush=True)
with open("e1_words_paid.txt", "w") as fh:
    fh.write("\n".join(keep) + "\n")
G = Game(FS)
for f in FS:
    t1 = time.time(); N = 4 ** f
    cols = []
    for w in keep:
        M = word_mat(G, w, f).tocoo()
        cols.append(sp.csr_matrix((M.data, (np.zeros(M.nnz, np.int64), M.row.astype(np.int64) * N + M.col)),
                                  shape=(1, N * N)))
    G.cache.clear()
    Wm = sp.vstack(cols).tocsr(); del cols
    GN = (Wm.conj() @ Wm.T).toarray()
    np.save("e1_gram_paid_f%d.npy" % f, GN)
    gf = GN / N
    ev = np.linalg.eigvalsh((gf + gf.conj().T) / 2)
    nz = ev[ev > 1e-9]
    print("f=%d N=%d: rank %d of %d; exact null %d; smallest nonzero eig %.4e; nonzero bottom: %s; %.1f s"
          % (f, N, nz.size, n, n - nz.size, nz.min(), " ".join("%.4f" % x for x in nz[:6]), time.time() - t1),
          flush=True)
print("done %.1f s" % (time.time() - t0))
