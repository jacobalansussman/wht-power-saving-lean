"""e1_analyze.py -- ex2-limits E1 analysis.  Reads e1_gram_<menu>_f<f>.npy (N*G_f) in the current folder.
For the column-uniform control: checks the tensor-power law N*G_f = (tr h_i^* h_j)^f via G_2 o G_3 = G_5, G_2 o G_2 = G_4.
For a menu: (i) rank and exact null space at each f; (ii) the relations exact at EVERY listed f >= fmin
(= null space of sum_f G_f, valid because each G_f is positive semidefinite); (iii) on the complement Q of those common
relations, the smallest eigenvalues of Q^* G_f Q at each f: asymptotic-only relations would show as eigenvalues that
decay towards 0 as f grows; (iv) entrywise convergence of G_f.
Usage: python3 e1_analyze.py ctl | <menu> <fmin> <f> <f> ...
"""
import sys, glob, numpy as np
menu = sys.argv[1]
def load(m, f):
    return np.load("e1_gram_%s_f%d.npy" % (m, f)) / 4.0 ** f
if menu == "ctl":
    G = {f: load("ctl", f) for f in (2, 3, 4, 5)}
    print("ctl: max |G2*G3 - G5| = %.2e ; max |G2*G2 - G4| = %.2e ; max |G5| offdiag = %.3f"
          % (np.abs(G[2] * G[3] - G[5]).max(), np.abs(G[2] * G[2] - G[4]).max(),
             np.abs(G[5] - np.diag(np.diag(G[5]))).max()))
    sys.exit()
fmin = int(sys.argv[2]); FS = [int(a) for a in sys.argv[3:]]
G = {f: load(menu, f) for f in FS}
n = G[FS[0]].shape[0]
for f in FS:
    ev = np.linalg.eigvalsh(G[f])
    nz = ev[ev > 1e-9]
    print("f=%d: rank %d of %d, per-f exact null %d, nonzero spectrum bottom: %s"
          % (f, nz.size, n, n - nz.size, " ".join("%.4f" % x for x in nz[:8])))
FC = [f for f in FS if f >= fmin]
S = sum(G[f] for f in FC)
w, U = np.linalg.eigh(S)
tol = 1e-9 * w.max()
Q = U[:, w > tol]; Z = U[:, w <= tol]
print("relations exact at every f in %s (fixed coefficients): %d ; complement dim %d" % (FC, Z.shape[1], Q.shape[1]))
for f in FS:
    ev = np.linalg.eigvalsh(Q.conj().T @ G[f] @ Q)
    print("  f=%d: on the complement: smallest eigenvalues %s ; # below 1e-9: %d"
          % (f, " ".join("%.5f" % x for x in ev[:8]), (ev < 1e-9).sum()))
for a, b in zip(FS[:-1], FS[1:]):
    print("  max |G_%d - G_%d| = %.3e" % (b, a, np.abs(G[b] - G[a]).max()))
# eigenvector of the smallest complement eigenvalue at the largest f: its words
fL = FS[-1]
ev, V = np.linalg.eigh(Q.conj().T @ G[fL] @ Q)
c = Q @ V[:, 0]
words = [l.strip() for l in open("e1_words_%s.txt" % menu)]
order = np.argsort(-np.abs(c))
print("smallest-eigenvalue combination at f=%d (eig %.5f), top terms:" % (fL, ev[0]))
for i in order[:10]:
    print("   %+.3f%+.3fi  %s" % (c[i].real, c[i].imag, words[i]))
print("  its norm^2 at each f:", " ".join("f%d:%.5f" % (f, (c.conj() @ G[f] @ c).real) for f in FS))
