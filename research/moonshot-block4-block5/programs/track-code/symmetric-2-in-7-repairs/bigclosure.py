"""bigclosure.py -- mw2-absorb4: full exact-structure test of a limit skeleton with a LARGE closure (sparse).
Closure over frames (bm.Game.closure, cap), then the sparse least-squares residual (scipy lsqr) of
[K]x_l = sum (column) beta at given content coefficients (the relaxation-exact point, and random perturbations of the
free coefficients that the relaxation leaves undetermined).  Usage: bigclosure.py lim.pkl cap npoints"""
import sys, pickle, time
from twist import *
import scipy.sparse as sp
from scipy.sparse.linalg import lsqr
d = pickle.load(open(sys.argv[1], 'rb'))
cap, npts = int(sys.argv[2]), int(sys.argv[3])
m, W, H = d['m'], d['W'], d['H']
sk = build_T(m, W, H, d['recipe'], d['skseed'], d.get('twist', {}))
G = sk.G
tok = []; gd = []
for i, T in enumerate(sk.gens):
    dd = {}
    for j, tt in enumerate(T):
        dd[tt] = len(tok) + 1; tok.append((i, j))
    gd.append(dd)
t0 = time.time()
cl = G.closure(gd, (0, G.KID), cap)
if cl is None:
    print(f"closure cap {cap} exceeded ({time.time()-t0:.0f}s)"); sys.exit()
rows, cols, ent = cl
print(f"closure rows {len(rows)} cols {len(cols)} entries {len(ent)} ({time.time()-t0:.0f}s)", flush=True)
E = np.array([(r, c, tok[co - 1][0], tok[co - 1][1]) for (r, c, co) in ent], dtype=np.int64)
sizes = [len(T) for T in sk.gens]; off = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
F = off[E[:, 2]] + E[:, 3]
N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]
al0 = np.asarray(d['alpha'], dtype=complex)
rng = np.random.default_rng(3)
for pt in range(npts):
    al = al0.copy()
    if pt > 0:
        al[free_idx] = al0[free_idx] * (1 + 0.3 * (rng.standard_normal(len(free_idx)) + 1j * rng.standard_normal(len(free_idx))))
    fl = np.concatenate(sk.coefs(al))
    A = sp.csr_matrix((fl[F], (E[:, 0], E[:, 1])), shape=(len(rows), len(cols)))
    t = np.zeros(len(rows), dtype=complex); t[0] = 1
    out = lsqr(A, t, atol=1e-14, btol=1e-14, iter_lim=20000)
    x = out[0]; r = A @ x - t
    print(f"point {pt}: residual^2 {np.sum(np.abs(r)**2):.4e} (lsqr stop {out[1]}, iters {out[2]}) ({time.time()-t0:.0f}s)", flush=True)
