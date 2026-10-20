# small11c.py -- 1+1 in 4 relay (hook-impossible) with pool letters: which two directions of the kernel's diagonal block
# are missed?  Prints the top residual eigenvectors as vectors over the 8 diagonal rows x (target rows K_x (x) K_x).
import numpy as np, sys
from fam import *
rng = np.random.default_rng(int(sys.argv[1]))
for tr in range(int(sys.argv[2])):
    sk = krelay(1, 1, 0, rng, 'struct', k=3)
    x, f = solve_bfgs(sk, rng, maxiter=600)
    c = x[:sk.n] + 1j * x[sk.n:]; G = sk.gens(c); W = sk.W
    Q = QB['D']
    S = np.vstack([np.hstack([Q @ g[l] for l in range(W)]) for g in G])
    U, s, Vh = np.linalg.svd(S, full_matrices=False); r = int((s > 1e-9 * s[0]).sum()); V = Vh[:r]
    T = np.zeros((8, 64 * W), complex); T[:, :64] = Q @ TG2['K']
    R = T - (T @ V.conj().T) @ V
    Gm = R @ R.conj().T; ev, vec = np.linalg.eigh(Gm)
    nrm = np.sum(np.abs(T) ** 2) / 8
    print('trial', tr, 'F', round(f, 4), 'kernel D-block Gram eigenvalues/row', np.round(ev[::-1][:3] / nrm, 3))
    for j in (-1, -2):
        v = vec[:, j]; v = v / v[np.argmax(np.abs(v))]
        print('   missed direction over x = 0..7:', np.round(v, 3))
