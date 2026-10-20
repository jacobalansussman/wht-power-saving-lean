import sys, numpy as np
menu = sys.argv[1]; fmin = int(sys.argv[2]); FS = [int(a) for a in sys.argv[3:]]
G = {f: np.load('my_G_paid_%s_f%d.npy' % (menu, f)) for f in FS}
G = {f: (g + g.conj().T) / 2 for f, g in G.items()}
S = sum(G[f] for f in FS if f >= fmin)
w, U = np.linalg.eigh(S); tol = 1e-9 * w.max(); Q = U[:, w > tol]
print('%s: relations exact at every f in %s: %d; complement dim %d' % (menu, [f for f in FS if f >= fmin], (w <= tol).sum(), Q.shape[1]))
for f in FS:
    ev = np.linalg.eigvalsh(Q.conj().T @ G[f] @ Q)
    print('  f=%d complement smallest: %s ; # < 1e-9: %d' % (f, ' '.join('%.5f' % x for x in ev[:6]), (ev < 1e-9).sum()))
