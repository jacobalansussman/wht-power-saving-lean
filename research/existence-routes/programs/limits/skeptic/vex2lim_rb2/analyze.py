"""analyze.py -- my analysis: relations exact at every f >= fmin (null space of sum of Grams), complement spectrum."""
import sys, numpy as np
menu = sys.argv[1]; fmin = int(sys.argv[2]); FS = [int(a) for a in sys.argv[3:]]
def load(f):
    c = np.load('my_cnt_%s_f%d.npy' % (menu, f)); return (c[0] - c[2] + 1j * (c[1] - c[3])) / 4.0 ** f
G = {f: load(f) for f in FS}
n = G[FS[0]].shape[0]
for f in FS:
    ev = np.linalg.eigvalsh(G[f]); print('f=%d rank %d/%d bottom nonzero %s' % (f, (ev > 1e-9).sum(), n, ' '.join('%.4f' % x for x in ev[ev > 1e-9][:5])))
S = sum(G[f] for f in FS if f >= fmin)
w, U = np.linalg.eigh(S); tol = 1e-9 * w.max()
Q = U[:, w > tol]
print('relations exact at every f in %s: %d; complement dim %d' % ([f for f in FS if f >= fmin], (w <= tol).sum(), Q.shape[1]))
for f in FS:
    ev = np.linalg.eigvalsh(Q.conj().T @ G[f] @ Q)
    print('  f=%d complement smallest: %s ; # < 1e-9: %d' % (f, ' '.join('%.5f' % x for x in ev[:6]), (ev < 1e-9).sum()))
